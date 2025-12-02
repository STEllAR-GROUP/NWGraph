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

    template <typename T>
    using local_pv_view = nw::graph::util::partitioned_vector_local_partition_view<T>;

    template <int idx, typename index_type, typename vertex_id, typename... Attributes>
    using local_adj_view =
      nw::graph::partitioned_index_adjacency_local_view<idx, index_type, vertex_id, Attributes...>;


    template <typename Graph, typename Real>
    using remote_results_type = std::vector<std::tuple<typename Graph::vertex_id_type, Real>>;

    template <typename Graph, typename Real, typename vertex_id_t = typename Graph::vertex_id_type>
    static remote_results_type<Graph, Real> do_page_rank_packet_4(
      std::size_t part_num, std::vector<std::tuple<vertex_id_t, vertex_id_t>>&& incoming_packet,
      hpx::partitioned_vector<Real> page_rank, hpx::partitioned_vector<vertex_id_t> degrees) {
      hpx::scoped_annotation _(__func__);

      // Construct the appropriate views for the remote action
      local_pv_view pr_view(page_rank, part_num);
      local_pv_view degrees_view(degrees, part_num);

      remote_results_type<Graph, Real> results;
      results.reserve(incoming_packet.size());
      for (auto&& [v_dest, v_asked] : incoming_packet) {
        if (degrees_view[v_asked] != 0) {
          Real outgoing_contribution = pr_view[v_asked] / degrees_view[v_asked];
          results.push_back(std::make_tuple(v_dest, outgoing_contribution));
        }
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

      template <typename ExPolicy, typename Graph, typename vertex_id_t = Graph::vertex_id_type>
      static Real sequential(ExPolicy&& policy, Graph G, size_t first_index, size_t last_index,
                             hpx::partitioned_vector<Real> page_rank_pv,
                             hpx::partitioned_vector<Real> accum_pv,
                             hpx::partitioned_vector<vertex_id_t> degrees_pv, Real base_score,
                             Real damping_factor, size_t batchsize) {
        // Convert to local views and call the main implementation
        std::size_t partnum = accum_pv.get_partition(first_index);
        local_pv_view page_rank_loc(page_rank_pv, partnum);
        local_pv_view degrees_loc(degrees_pv, partnum);
        local_pv_view accum_loc(accum_pv, partnum);
        local_adj_view G_loc(G, partnum);
        return seq_impl(policy, G_loc, page_rank_loc, accum_loc, degrees_loc, base_score,
                        damping_factor, batchsize);
      }

      template <typename ExPolicy, typename LocGraph,
                typename vertex_id_t = typename LocGraph::vertex_id_type>
      static Real seq_impl(ExPolicy&& policy, LocGraph G_loc, local_pv_view<Real> page_rank_loc,
                           local_pv_view<Real> accum_loc, local_pv_view<vertex_id_t> degrees_loc,
                           Real base_score, Real damping_factor, size_t batchsize) {
        hpx::scoped_annotation _(__func__);

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

        using Graph = LocGraph::graph_type;
        using loc_iter_t = LocGraph::iterator;

        auto get_part_locality = [](auto& G_loc, std::size_t part_num) -> hpx::id_type
        {
          // TODO: Too intrusive, fix
          hpx::id_type part_id = G_loc.parent().get_indices().partitions()[part_num].get_id();
          return hpx::naming::get_locality_from_id(part_id);
        };

        auto send_remote_action = [&](std::size_t part_num, auto&& targets) -> hpx::future<void>
        {
          // hpx::scoped_annotation __("PR4:send_remote_action");
          using action_t = page_rank_action_4<Graph, Real>;
          // Get locality of target partition
          hpx::id_type locality_id = get_part_locality(G_loc, part_num);

          auto f1 =
            hpx::async<action_t>(locality_id, part_num, std::move(targets),
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
        std::map<part_num_t, std::vector<std::tuple<vertex_id_t, vertex_id_t>>> outgoing_packets;

        std::vector<hpx::future<void>> remote_results;
        std::vector<std::tuple<vertex_id_t, vertex_id_t>> v_targets;

        {
          hpx::scoped_annotation __("gather contributions");
          for (loc_iter_t v_it = G_loc.begin(); v_it != G_loc.end(); ++v_it) {
            // for (auto&& edge_rng : G_loc) {
            auto edge_rng = *v_it;
            vertex_id_t u = v_it.index();

            for (auto&& edge : edge_rng) {

              vertex_id_t v = target(G_loc, edge);

              if (page_rank_loc.is_local_index(v) && degrees_loc[v] != 0) {
                std::atomic_ref(accum_loc[u]) += page_rank_loc[v] / degrees_loc[v];
              }
              else {
                auto part_num = vertex_partition_num(G_loc.parent(), v);
                outgoing_packets[part_num].push_back({u, v});
              }
            }

            for (auto& [part_num, targets] : outgoing_packets) {
              if (targets.size() > batchsize) {
                remote_results.emplace_back(send_remote_action(part_num, std::move(targets)));
                targets = {};
              }
            }
          }

          // send any remaining requests
          for (auto& [part_num, targets] : outgoing_packets) {
            if (!targets.empty()) {
              remote_results.emplace_back(send_remote_action(part_num, std::move(targets)));
            }
          }
        }

        // wait for all remote operations to finish
        hpx::wait_all(remote_results);

        // (local) accumulated result is now fully computed
        // update the local page rank
        Real local_error = 0.0;
        {
          hpx::scoped_annotation __("update page rank");
          for (loc_iter_t v_it = G_loc.begin(); v_it != G_loc.end(); ++v_it) {
            vertex_id_t u = v_it.index();

            Real old_rank = page_rank_loc[u];
            Real new_rank = base_score + damping_factor * accum_loc[u];

            page_rank_loc[u] = new_rank;
            local_error += fabs(new_rank - old_rank);
            accum_loc[u] = 0.0;
          }
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
    hpx::fill(hpx::execution::par_unseq, page_rank.begin(), page_rank.end(), init_score);

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


      std::vector<hpx::future<Real>> errors = partitioned_algorithm<detail::page_rank_4<Real>>(
        hpx::execution::seq, G, hpx::ref(page_rank), hpx::ref(accummulating_contributions),
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
