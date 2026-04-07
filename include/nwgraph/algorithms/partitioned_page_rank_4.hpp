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

#include "nwgraph/algorithms/partitioned_algorithm.hpp"


#include <algorithm>
#include <chrono>
#include <cstddef>
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


    template <typename Graph, typename Real>
    using remote_results_type =
      std::vector<std::tuple<vertex_id_t<std::remove_reference_t<Graph>>, Real>>;

    template <typename Graph, typename Partition, typename Real>
    static remote_results_type<Graph, Real> do_page_rank_packet_4(
      Partition partition,
      std::vector<std::tuple<vertex_id_t<std::remove_reference_t<Graph>>,
                             vertex_id_t<std::remove_reference_t<Graph>>>>&& incoming_packet,
      hpx::partitioned_vector<Real> page_rank,
      hpx::partitioned_vector<vertex_id_t<std::remove_reference_t<Graph>>> degrees) {

      hpx::scoped_annotation ann_pr_handle_packet("PR_handle_packet");

      decltype(auto) pr_data = detail::partition_data(page_rank, partition);
      decltype(auto) degrees_data = detail::partition_data(degrees, partition);

      remote_results_type<Graph, Real> results;
      results.reserve(incoming_packet.size());
      for (auto&& [v_dest, v_asked] : incoming_packet) {
        if (degrees_data[v_asked] != 0) {
          Real outgoing_contribution = pr_data[v_asked] / degrees_data[v_asked];
          results.push_back(std::make_tuple(v_dest, outgoing_contribution));
        }
      }
      return results;
    }

    template <typename Graph, typename Partition, typename Real>
    struct page_rank_action_4
      : hpx::actions::action<decltype(&do_page_rank_packet_4<Graph, Partition, Real>),
                             &do_page_rank_packet_4<Graph, Partition, Real>,
                             page_rank_action_4<Graph, Partition, Real>> {
    };

    ////////////////////////////////////////////////////////////////////////////
    template <typename Real>
    struct page_rank_4 : hpx::parallel::detail::algorithm<page_rank_4<Real>, Real> {
      static constexpr bool partition_aware = true;

      // page rank driver for one of the partitions
      constexpr page_rank_4() noexcept
        : hpx::parallel::detail::algorithm<page_rank_4, Real>("page_rank_4") {}

      template <typename ExPolicy, typename Graph>
      static Real sequential(ExPolicy&& policy, Graph G, partition_t<std::remove_reference_t<Graph>> this_partition,
                             hpx::partitioned_vector<Real> page_rank_pv,
                             hpx::partitioned_vector<Real> accum_pv,
                             hpx::partitioned_vector<vertex_id_t<std::remove_reference_t<Graph>>> degrees_pv, Real base_score,
                             Real damping_factor, size_t batchsize) {
        hpx::scoped_annotation ann_pr_impl("PR_impl");

        decltype(auto) G_data = detail::partition_data(G, this_partition);
        decltype(auto) page_rank_data = detail::partition_data(page_rank_pv, this_partition);
        decltype(auto) degrees_data = detail::partition_data(degrees_pv, this_partition);
        decltype(auto) accum_data = detail::partition_data(accum_pv, this_partition);

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
        auto first_index = static_cast<vertex_id_type>(partition_first_index(this_partition));
        auto last_index = static_cast<vertex_id_type>(partition_last_index(this_partition));

        auto send_remote_action = [&](partition_type const& target_partition,
                                      auto&& targets) -> hpx::future<void>
        {
          using action_t = page_rank_action_4<graph_type, partition_type, Real>;

          auto f1 = hpx::async<action_t>(
            partition_locality(target_partition), target_partition, std::move(targets),
            hpx::ref(page_rank_pv), hpx::ref(degrees_pv));
          return f1.then(
            [&](auto&& f)
            {
              for (auto&& [v, incoming_val] : f.get()) {
                std::atomic_ref(accum_data[v]) += incoming_val;
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
          first_index,
          last_index,
          [&](vertex_id_type u)
                      {
                        auto* thd_outgoing_packets = &(outgoing_packets.get());
                        auto edge_rng = detail::iterable_row(G_data[u]);

                        for (auto&& edge : edge_rng) {
                          vertex_id_type v = detail::edge_target(G_data, edge);

                          if (is_local(this_partition, v)) {
                            if (degrees_data[v] != 0)
                              std::atomic_ref(accum_data[u]) += page_rank_data[v] / degrees_data[v];
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
                                      first_index,
                                      last_index,
                                      hpx::experimental::reduction_plus(error),
                                      [&](vertex_id_type u, auto& local_error)
                                      {
                                        Real old_rank = page_rank_data[u];
                                        Real new_rank = base_score + damping_factor * accum_data[u];

                                        page_rank_data[u] = new_rank;
                                        local_error += fabs(new_rank - old_rank);
                                        accum_data[u] = 0.0;
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
