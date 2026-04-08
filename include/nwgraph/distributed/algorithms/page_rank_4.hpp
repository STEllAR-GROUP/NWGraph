/**
 * @file partitioned_page_rank_4.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * @authors
 *   Hartmut Kaiser
 *   Panagiotis (Panos) Syskakis
 */

#ifndef NW_GRAPH_PARTITIONED_PAGE_RANK_4_HPP
#define NW_GRAPH_PARTITIONED_PAGE_RANK_4_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/distributed/algorithms/algorithm.hpp"
#include "nwgraph/distributed/copartitioned_vector.hpp"


#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cmath>
#include <map>
#include <vector>

#include <hpx/async_combinators/wait_all.hpp>
#include <hpx/executors/execution_policy.hpp>
#include <hpx/include/partitioned_vector_predef.hpp>
#include <hpx/parallel/segmented_algorithms/detail/dispatch.hpp>
#include <hpx/parallel/util/detail/algorithm_result.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace nw::graph {

  namespace detail {

    namespace page_rank_4_impl {


    template <typename Graph, typename Real>
    using remote_results_type =
      std::vector<std::tuple<vertex_id_t<std::remove_reference_t<Graph>>, Real>>;

    template <typename Graph, typename Real>
    static remote_results_type<Graph, Real> do_page_rank_packet_4(
      std::vector<std::tuple<vertex_id_t<std::remove_reference_t<Graph>>,
                             vertex_id_t<std::remove_reference_t<Graph>>>>&& incoming_packet,
      hpx::partitioned_vector<Real> page_rank,
      hpx::partitioned_vector<vertex_id_t<std::remove_reference_t<Graph>>> degrees) {

      hpx::scoped_annotation ann_pr_handle_packet("PR_handle_packet");

      remote_results_type<Graph, Real> results;
      results.reserve(incoming_packet.size());
      for (auto&& [v_dest, v_asked] : incoming_packet) {
        auto degree = *degrees.get_local_iterator(static_cast<std::size_t>(v_asked)).local();
        if (degree != 0) {
          auto rank = *page_rank.get_local_iterator(static_cast<std::size_t>(v_asked)).local();
          Real outgoing_contribution = rank / degree;
          results.push_back(std::make_tuple(v_dest, outgoing_contribution));
        }
      }
      return results;
    }

    template <typename Graph, typename Real>
    struct page_rank_action_4
      : hpx::actions::action<decltype(&do_page_rank_packet_4<Graph, Real>),
                             &do_page_rank_packet_4<Graph, Real>,
                             page_rank_action_4<Graph, Real>> {
    };

    } // namespace page_rank_4_impl

    ////////////////////////////////////////////////////////////////////////////
    template <typename Real>
    struct page_rank_4 : hpx::parallel::detail::algorithm<page_rank_4<Real>, Real> {
      // page rank driver for one of the partitions
      constexpr page_rank_4() noexcept
        : hpx::parallel::detail::algorithm<page_rank_4, Real>("page_rank_4") {}

      template <typename ExPolicy, typename Graph>
      static Real sequential(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index,
                             hpx::partitioned_vector<Real> page_rank_pv,
                             hpx::partitioned_vector<Real> accum_pv,
                             hpx::partitioned_vector<vertex_id_t<std::remove_reference_t<Graph>>> degrees_pv, Real base_score,
                             Real damping_factor, size_t batchsize) {
        hpx::scoped_annotation ann_pr_impl("PR_impl");

        // Compute one Page-Rank iteration
        // For each local node, do:
        // 1. Get contribution from neighbors. For each vertex, compute the outgoing contribution as
        // page_rank[i] / degrees[i] 1.1 For remote neighbors, send action to each locality that has
        // at least one neighbor, and compute the contribution of that locality on this node. Send
        // the result back to this locality. 1.2 For local neighbors, update contribution directly
        // 2. Update page rank for each vertex using the formula: page_rank[i] = base_score +
        // damping_factor * z where z is the sum of contributions from all neighbors
        // 3. Return the local accumulated error (sum of absolute difference between old and new
        // page ranks)

        // Lambda to request contributions from remote localities, and update the local page rank
        // once the contributions are received

        using graph_type = std::remove_reference_t<Graph>;
        using vertex_id_type = vertex_id_t<graph_type>;
        using partition_type = partition_t<graph_type>;
        auto first_vertex = static_cast<vertex_id_type>(first_index);
        auto last_vertex = static_cast<vertex_id_type>(last_index);

        auto rank_first = page_rank_pv.get_local_iterator(first_index).local();
        auto accum_first = accum_pv.get_local_iterator(first_index).local();
        auto degree_first = degrees_pv.get_local_iterator(first_index).local();

        auto send_remote_action = [&](partition_type const& target_partition,
                                      auto&& targets) -> hpx::future<void>
        {
          using action_t = detail::page_rank_4_impl::page_rank_action_4<graph_type, Real>;

          auto f1 = hpx::async<action_t>(
            static_cast<hpx::id_type>(target_partition), std::move(targets),
            hpx::ref(page_rank_pv), hpx::ref(degrees_pv));
          return f1.then(
            [&](auto&& f)
            {
              for (auto&& [v, incoming_val] : f.get()) {
                auto& accum = *(accum_pv.get_local_iterator(static_cast<std::size_t>(v)).local());
                std::atomic_ref(accum) += incoming_val;
              }
            });
        };

        using packet_t = std::vector<std::tuple<vertex_id_type, vertex_id_type>>;
        safe_object<std::map<partition_type, packet_t>> outgoing_packets;

        using results_t = std::vector<hpx::future<void>>;
        safe_object<results_t> remote_results;


        {
          hpx::scoped_annotation ann_pr_for_each("PR_for_each");

        hpx::experimental::for_loop(
          hpx::execution::par,
          first_vertex,
          last_vertex,
          [&](vertex_id_type u)
                      {
                        auto* thd_outgoing_packets = &(outgoing_packets.get());
                        auto edge_rng = G[u];
                        auto u_offset = static_cast<std::size_t>(u - first_vertex);

                        for (auto&& edge : edge_rng) {
                          vertex_id_type v = static_cast<vertex_id_type>(target(G, edge));

                          if (first_vertex <= v && v < last_vertex) {
                            auto v_offset = static_cast<std::size_t>(v - first_vertex);
                            auto& degree = *(degree_first + v_offset);
                            if (degree != 0) {
                              auto& rank = *(rank_first + v_offset);
                              auto& accum = *(accum_first + u_offset);
                              std::atomic_ref(accum) += rank / degree;
                            }
                          }
                          else {
                            auto target_partition = vertex_partition(G, v);
                            (*thd_outgoing_packets)[target_partition].push_back({u, v});
                          }
                        }

                        for (auto& [target_partition, targets] : *thd_outgoing_packets) {
                          // If the hpx thread ever suspended, it could have migrated to another
                          // thread, in which case it is no longer safe to use the thread-local
                          // reference
                          if (thd_outgoing_packets != &(outgoing_packets.get()))
                            break;
                          if (targets.size() > batchsize) {
                            auto tmp = std::move(targets);
                            targets = {};
                            auto f = send_remote_action(target_partition, std::move(tmp));
                            remote_results.get().emplace_back(std::move(f));
                          }
                        }
                      });

        }


        // send any remaining requests
        outgoing_packets.reduce(
          [&](auto&& thd_outgoing_packets)
          {
            for (auto& [target_partition, targets] : thd_outgoing_packets) {
              if (!targets.empty()) {
                remote_results.get().emplace_back(
                  send_remote_action(target_partition, std::move(targets)));
              }
            }
          });


        // wait for all remote operations to finish
        {
          hpx::scoped_annotation ann_pr_wait_all("PR_wait_all"); 
          remote_results.reduce([](auto&& thd_remote_results) { hpx::wait_all(thd_remote_results); });
        }

        // (local) accumulated result is now fully computed
        // update the local page rank
        Real error = 0.0;

        {
          hpx::scoped_annotation ann_pr_update_pr("PR_update_pr");

          hpx::experimental::for_loop(
                                      hpx::execution::par,
                                      first_vertex,
                                      last_vertex,
                                      hpx::experimental::reduction_plus(error),
                                      [&](vertex_id_type u, auto& local_error)
                                      {
                                        auto offset = static_cast<std::size_t>(u - first_vertex);
                                        auto& rank = *(rank_first + offset);
                                        auto& accum = *(accum_first + offset);
                                        Real old_rank = rank;
                                        Real new_rank = base_score + damping_factor * accum;

                                        rank = new_rank;
                                        local_error += fabs(new_rank - old_rank);
                                        accum = 0.0;
                                      });
        }

        return error;
      }


      template <typename ExPolicy, typename Graph, typename IterB, typename IterE>
      static size_t parallel(ExPolicy&& policy, Graph G, IterB first, IterE last) {
        return 0;
      }
    };


    /// \endcond
  } // namespace detail

  template <partitioned_algorithm_graph Graph, typename Real = double>
  void partitioned_page_rank_4(Graph& G,
                               hpx::partitioned_vector<vertex_id_t<std::remove_reference_t<Graph>>>& degrees,
                               hpx::partitioned_vector<Real>& page_rank,
                               const Real damping_factor = 0.85, const Real threshold = 1.e-4,
                               const size_t max_iters = std::numeric_limits<unsigned int>::max(),
                               size_t batchsize = 1000) {

    hpx::scoped_annotation ann_pr("PR");

    const Real init_score = 1.0 / G.size();
    const Real base_score = (1.0 - damping_factor) / G.size();

    // Initialize page ranks
    hpx::fill(hpx::execution::par_unseq, page_rank.begin(), page_rank.end(), init_score);

    auto accumulating_contributions = make_copartitioned_vector<Real>(G, 0.0);
    accumulating_contributions.register_as("accumulating_contributions");

    for (size_t iter = 0; iter < max_iters; ++iter) {

      hpx::scoped_annotation ann_pr_iteration("PR_iteration");

      std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
      std::cout << "----- Iteration " << iter << " ----- " << std::endl;

      std::vector<hpx::future<Real>> errors = partitioned_algorithm<detail::page_rank_4<Real>>(
        hpx::execution::seq, G, hpx::ref(page_rank), hpx::ref(accumulating_contributions),
        hpx::ref(degrees), base_score, damping_factor, batchsize);

      // std::transform_reduce is undefined behavior, since it modifies the input futures
      Real error = 0.0;
      for (auto& f : errors) {
        error += f.get();
      }

      std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed_seconds = end_time - start_time;
      std::cout << "Iteration " << iter << " completed in " << elapsed_seconds.count()
                << " seconds with error " << error << std::endl;
      if (error < threshold)
        break;
    }
  }
} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_PAGE_RANK_4_HPP
