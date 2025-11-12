/**
 * @file partitioned_page_rank_1.hpp
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

#ifndef NW_GRAPH_PARTITIONED_PAGE_RANK_1_HPP
#define NW_GRAPH_PARTITIONED_PAGE_RANK_1_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include <nwgraph/util/partitioned_vector_local_partition_view.hpp>
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
    using remote_results_type = std::vector<std::tuple<typename Graph::vertex_id_type, Real>>;

    template <typename Graph, typename Real>
    static remote_results_type<Graph, Real> do_page_rank_packet_1(
      std::vector<std::tuple<typename Graph::vertex_id_type, typename Graph::vertex_id_type>>&&
        incoming_packet,
      hpx::partitioned_vector<Real> page_rank,
      hpx::partitioned_vector<typename Graph::vertex_id_type> degrees) {

      remote_results_type<Graph, Real> results;
      for (auto&& [v_dest, v_asked] : incoming_packet) {
        auto deg_iter = degrees.get_local_iterator(v_asked).local();
        auto pr_iter = page_rank.get_local_iterator(v_asked).local();
        Real outgoing_contribution = *pr_iter / *deg_iter;
        results.push_back(std::make_tuple(v_dest, outgoing_contribution));
      }
      return results;
    }

    template <typename Graph, typename Real>
    struct page_rank_action_1
      : hpx::actions::action<decltype(&do_page_rank_packet_1<Graph, Real>),
                             &do_page_rank_packet_1<Graph, Real>, page_rank_action_1<Graph, Real>> {
    };

    ////////////////////////////////////////////////////////////////////////////
    template <typename Real>
    struct page_rank_1 : hpx::parallel::detail::algorithm<page_rank_1<Real>, Real> {

      // page rank driver for one of the partitions
      constexpr page_rank_1() noexcept
        : hpx::parallel::detail::algorithm<page_rank_1, Real>("page_rank_1") {}

      template <typename T>
      using local_partition_view =
        typename nw::graph::util::partitioned_vector_local_partition_view<T, std::vector<T>>;

      template <typename ExPolicy, typename Graph>
      static Real sequential(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index, hpx::partitioned_vector<Real> page_rank_pv,
                             hpx::partitioned_vector<Real> accum_pv,
                             hpx::partitioned_vector<typename Graph::vertex_id_type> degrees_pv,
                             Real base_score, Real damping_factor, size_t batchsize) {

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

        using vertex_id_type = typename Graph::vertex_id_type;

        local_partition_view<Real> page_rank_loc(page_rank_pv, page_rank_pv.get_partition(first_index));
        local_partition_view<vertex_id_type> degrees_loc(degrees_pv, degrees_pv.get_partition(first_index));
        local_partition_view<Real> accum_loc(accum_pv, accum_pv.get_partition(first_index));

        auto send_request_packet =
          [&page_rank_pv, &degrees_pv, &accum_pv](
            auto&& id, auto&& targets) -> hpx::future<remote_results_type<Graph, Real>>
        {
          using action_t = page_rank_action_1<Graph, Real>;
          return hpx::async<action_t>(id, std::move(targets), hpx::ref(page_rank_pv),
                                      hpx::ref(degrees_pv));
        };

        auto handle_response = [&accum_pv](auto&& f) -> void
        {
          using remote_result_t = remote_results_type<Graph, Real>;
          remote_result_t incoming_packet = f.get();
          for (auto&& [v, incoming_val] : incoming_packet) {
            auto acc_iter = accum_pv.get_local_iterator(v).local();
            std::atomic_ref(*acc_iter) += incoming_val;
          }
        };


        hpx::id_type this_locality_id = hpx::find_here();

        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;
        auto G_sub = std::ranges::subrange(first, last);


        std::map<hpx::id_type, std::vector<std::tuple<vertex_id_type, vertex_id_type>>>
          outgoing_packets;

        std::vector<hpx::future<void>> remote_results;
        std::vector<std::tuple<vertex_id_type, vertex_id_type>> v_targets;

        // for each v in G do
        // for (auto edge_rng : G_sub) {
        for (auto v_it = first; v_it != last; ++v_it) {
          auto edge_rng = *v_it;
          for (auto edge : edge_rng) {

            // vertex_id_type u = source(G, edge);
            //     TODO: getting "u" looks ugly, we could keep the src vertex in the edge range or
            //     the edge, or create a nicer getter function for it
            vertex_id_type u = v_it.index().get_global_index();
            vertex_id_type v = target(G, edge);

            auto loc_id = vertex_locality(G, v);
            if (this_locality_id == loc_id) {
              // handle things locally
              accum_loc[u] += page_rank_loc[v] / degrees_loc[v];
            }
            else {
              outgoing_packets[loc_id].push_back({u, v});
            }
          }

          for (auto&& [id, targets] : outgoing_packets) {
            if (targets.size() < batchsize)
              continue;

            remote_results.emplace_back(
              send_request_packet(id, std::move(targets)).then(handle_response));
            targets = {};
          }
        }

        // send any remaining requests
        for (auto&& [id, targets] : outgoing_packets) {
          if (targets.empty())
            continue;
          remote_results.emplace_back(
            send_request_packet(id, std::move(targets)).then(handle_response));
        }

        // wait for all remote operations to finish
        hpx::wait_all(remote_results);

        // (local) accumulated result is now fully computed
        // update the local page rank
        Real local_error = 0.0;


        // vertex_id_type i = first_index;

        for (auto v_it = first; v_it != last; ++v_it) {
          vertex_id_type u = v_it.index().get_global_index();

          Real old_rank = page_rank_loc[u];
          Real new_rank = base_score + damping_factor * accum_loc[u];
          
          page_rank_loc[u] = new_rank;
          local_error += fabs(new_rank - old_rank);
          accum_loc[u] = 0.0;
        }

        return local_error;
      }

      template <typename ExPolicy, typename Graph, typename IterB, typename IterE>
      static size_t parallel(ExPolicy&& policy, Graph G, IterB first, IterE last) {
        return 0;
      }
    };


    /// \endcond
  } // namespace detail

  template <adjacency_list_graph Graph, typename Real = double>
  void partitioned_page_rank_1(Graph& G,
                               hpx::partitioned_vector<typename Graph::vertex_id_type>& degrees,
                               hpx::partitioned_vector<Real>& page_rank,
                               const Real damping_factor = 0.85, const Real threshold = 1.e-4,
                               const size_t max_iters = std::numeric_limits<unsigned int>::max(),
                               size_t batchsize = 1000) {

    auto sizes = G.indices_.get_partition_sizes();

    const Real init_score = 1.0 / G.size();
    const Real base_score = (1.0 - damping_factor) / G.size();

    // Initialize page ranks
    hpx::fill(hpx::execution::seq, page_rank.begin(), page_rank.end(), init_score);

    hpx::partitioned_vector<Real> accummulating_contributions(
      G.indices_.size(), 0.0,
      hpx::explicit_container_layout(sizes, G.indices_.get_partition_localities()));
    accummulating_contributions.register_as("accummulating_contributions");

    for (size_t iter = 0; iter < max_iters; ++iter) {

      std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
      std::cout << "----- Iteration " << iter << " ----- " << std::endl;

      // for (auto i = 0; i < G.size(); ++i) {
      //   std::cout << "Node " << i << " : " << page_rank[i] << std::endl;
      // }


      auto errors = partitioned_algorithm<detail::page_rank_1<Real>>(
        hpx::execution::seq, G, hpx::ref(page_rank), hpx::ref(accummulating_contributions),
        hpx::ref(degrees), base_score, damping_factor, batchsize);

      auto error = std::transform_reduce(
        errors.begin(), errors.end(), Real(0.0), [](Real count, Real curr) { return count + curr; },
        [](auto&& f) { return f.get(); });

      std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed_seconds = end_time - start_time;
      std::cout << "Iteration " << iter << " completed in " << elapsed_seconds.count()
                << " seconds with error " << error << std::endl;
      if (error < threshold)
        break;
    }
  }
} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_PAGE_RANK_1_HPP
