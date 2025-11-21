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

#include <nwgraph/containers/partitioned_compressed_local_view.hpp>
#include <nwgraph/partitioned_adjacency_local_view.hpp>
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

    // TODO: Where do these go?
    template <typename T>
    using local_partition_view =
      typename nw::graph::util::partitioned_vector_local_partition_view<T>;

    template <typename T>
    static auto to_local_view(hpx::partitioned_vector<T>& pv, std::size_t partnum) {
      return local_partition_view<T>(pv, partnum);
    }

    template <int idx, typename... Attributes>
    using local_adjacency_view = nw::graph::partitioned_adjacency_local_view<idx, Attributes...>;

    template <int idx, typename... Attributes>
    static auto to_local_view(nw::graph::partitioned_adjacency<idx, Attributes...>& G,
                              std::size_t partnum) {
      return local_adjacency_view<idx, Attributes...>(G, partnum);
    }


    template <typename Graph, typename Real>
    using remote_results_type = std::vector<std::tuple<typename Graph::vertex_id_type, Real>>;

    template <typename Graph, typename Real>
    static remote_results_type<Graph, Real> do_page_rank_packet_4(
      std::size_t part_num,
      std::vector<std::tuple<typename Graph::vertex_id_type, typename Graph::vertex_id_type>> &&
        incoming_packet,
      hpx::partitioned_vector<Real> page_rank,
      hpx::partitioned_vector<typename Graph::vertex_id_type> degrees) {

      // Construct the appropriate views for the remote action
      auto pr_view = to_local_view(page_rank, part_num);
      auto degrees_view = to_local_view(degrees, part_num);

      remote_results_type<Graph, Real> results;
      for (auto&& [v_dest, v_asked] : incoming_packet) {
        Real outgoing_contribution = pr_view[v_asked] / degrees_view[v_asked];
        results.push_back(std::make_tuple(v_dest, outgoing_contribution));
      }
      return results;
    }

    template <typename Graph, typename Real>
    struct page_rank_action_4
      : hpx::actions::action<decltype(&do_page_rank_packet_4<Graph, Real>),
                             &do_page_rank_packet_4<Graph, Real>, page_rank_action_4<Graph, Real>> {
    };

    ////////////////////////////////////////////////////////////////////////////
    template <typename Real>
    struct page_rank_4 : hpx::parallel::detail::algorithm<page_rank_4<Real>, Real> {

      // page rank driver for one of the partitions
      constexpr page_rank_4() noexcept
        : hpx::parallel::detail::algorithm<page_rank_4, Real>("page_rank_4") {}

      template <typename ExPolicy, typename Graph>
      static Real sequential(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index, hpx::partitioned_vector<Real> page_rank_pv,
                             hpx::partitioned_vector<Real> accum_pv,
                             hpx::partitioned_vector<typename Graph::vertex_id_type> degrees_pv,
                             Real base_score, Real damping_factor, size_t batchsize) {
        // Convert to local views and call the main implementation
        std::size_t partnum = accum_pv.get_partition(first_index);
        auto page_rank_loc = to_local_view(page_rank_pv, partnum);
        auto degrees_loc = to_local_view(degrees_pv, partnum);
        auto accum_loc = to_local_view(accum_pv, partnum);
        auto G_loc = to_local_view(G, partnum);
        return seq_impl(policy, G_loc, page_rank_loc, accum_loc, degrees_loc, base_score,
                        damping_factor, batchsize);
      }

      template <typename ExPolicy, int G_idx, typename... G_Attributes>
      static Real seq_impl(
        ExPolicy&& policy, local_adjacency_view<G_idx, G_Attributes...> G_loc,
        local_partition_view<Real> page_rank_loc, local_partition_view<Real> accum_loc,
        local_partition_view<typename local_adjacency_view<G_idx, G_Attributes...>::vertex_id_type>
          degrees_loc,
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

        using vertex_id_type = local_adjacency_view<G_idx, G_Attributes...>::vertex_id_type;
        using Graph = local_adjacency_view<G_idx, G_Attributes...>::graph_type;
        using loc_iter_t = local_adjacency_view<G_idx, G_Attributes...>::iterator;

        auto send_remote_action = [&](std::size_t part_num, auto&& targets) -> hpx::future<void>
        {
          using action_t = page_rank_action_4<Graph, Real>;
          // Get locality of target partition
          // TODO: Find a nicer way to do this (make locality accessible from X_view, perhaps
          // also make it support hpx::colocated)
          hpx::id_type part_id = page_rank_loc.parent().partitions()[part_num].get_id();
          hpx::id_type locality_id = hpx::naming::get_locality_from_id(part_id);
          auto f1 =
            hpx::async<action_t>(locality_id, part_num, std::move(targets), // TODO move
                                 hpx::ref(page_rank_loc.parent()), hpx::ref(degrees_loc.parent()));
          return f1.then(
            [&](auto&& f)
            {
              for (auto&& [v, incoming_val] : f.get()) {
                std::atomic_ref(accum_loc[v]) += incoming_val;
              }
            });
        };

        using part_num_t = std::size_t;
        std::map<part_num_t, std::vector<std::tuple<vertex_id_type, vertex_id_type>>>
          outgoing_packets;

        std::vector<hpx::future<void>> remote_results;
        std::vector<std::tuple<vertex_id_type, vertex_id_type>> v_targets;

        for (loc_iter_t v_it = G_loc.begin(); v_it != G_loc.end(); ++v_it) {
          // for (auto&& edge_rng : G_loc) {
          auto edge_rng = *v_it;
          vertex_id_type u = v_it.index();

          for (auto&& edge : edge_rng) {

            // vertex_id_type v = target(G_loc, edge);
            vertex_id_type v = std::get<0>(edge); // TODO: fullfill CPO concepts

            if (page_rank_loc.is_local_index(v)) {
              if (degrees_loc[v] == 0)
                continue;
              accum_loc[u] += page_rank_loc[v] / degrees_loc[v];
            }
            else {
              auto part_num = vertex_partition_num(G_loc.parent(), v);
              outgoing_packets[part_num].push_back({u, v});
            }
          }

          for (auto&& [part_num, targets] : outgoing_packets) {
            if (targets.size() > batchsize) {
              remote_results.emplace_back(send_remote_action(part_num, std::move(targets)));
              targets = {};
            }
          }
        }

        // send any remaining requests
        for (auto&& [part_num, targets] : outgoing_packets) {
          if (!targets.empty()) {
            remote_results.emplace_back(send_remote_action(part_num, std::move(targets)));
          }
        }

        // wait for all remote operations to finish
        hpx::wait_all(remote_results);

        // (local) accumulated result is now fully computed
        // update the local page rank
        Real local_error = 0.0;

        for (loc_iter_t v_it = G_loc.begin(); v_it != G_loc.end(); ++v_it) {
          vertex_id_type u = v_it.index();

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
  void partitioned_page_rank_4(Graph& G,
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


      auto errors = partitioned_algorithm<detail::page_rank_4<Real>>(
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

#endif // NW_GRAPH_PARTITIONED_PAGE_RANK_4_HPP
