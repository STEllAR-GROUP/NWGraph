/**
 * @file partitioned_page_rank_b.hpp
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

// b stands for backwards, because I misinterpreted the direction of the arrows in the graph :(

#ifndef NW_GRAPH_PARTITIONED_PAGE_RANK_B_HPP
#define NW_GRAPH_PARTITIONED_PAGE_RANK_B_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/distributed/algorithms/algorithm.hpp"

#include <algorithm>
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

    ////////////////////////////////////////////////////////////////////////////
    // handle incoming packets
    template <typename Graph, typename Real>
    static void page_rank_packet_b_do_remote_packet(
      Graph G, std::vector<std::tuple<typename Graph::vertex_id_type, Real>> const& incoming_packet,
      hpx::partitioned_vector<Real> accumulated_contributions) {

      for (auto&& [v, incoming_val] : incoming_packet) {
        accumulated_contributions[v] = incoming_val + accumulated_contributions[v];
      }
    }

    template <typename Graph, typename Real>
    struct page_rank_action_b
      : hpx::actions::action<decltype(&page_rank_packet_b_do_remote_packet<Graph, Real>),
                             &page_rank_packet_b_do_remote_packet<Graph, Real>,
                             page_rank_action_b<Graph, Real>> {
    };

    ////////////////////////////////////////////////////////////////////////////
    struct page_rank_0 : hpx::parallel::detail::algorithm<page_rank_0, size_t> {

      // page rank driver for one of the partitions
      constexpr page_rank_0() noexcept
        : hpx::parallel::detail::algorithm<page_rank_0, size_t>("page_rank_0") {}

      template <typename ExPolicy, typename Graph, typename Real>
      static size_t sequential(ExPolicy&& policy, Graph G, const size_t first_index,
                               const size_t last_index,
                               hpx::partitioned_vector<typename Real> page_rank,
                               hpx::partitioned_vector<typename Real> accumulated_contributions,
                               hpx::partitioned_vector<typename Graph::vertex_id_type> degrees) {


        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;

        std::uint32_t this_locality_id = hpx::get_locality_id();

        using vertex_id_type = typename Graph::vertex_id_type;


        std::map<hpx::id_type, std::vector<std::tuple<vertex_id_type, Real>>> outgoing_packets;
        std::vector<hpx::future<void>> remote_ops;
        std::vector<vertex_id_type> v_targets;

        // for each v in G do
        for (auto v_it = first; v_it != last; ++v_it) {

          v_targets.resize(0);

          size_t v_idx = first_index + std::distance(first, v_it);

          if (degrees[v_idx] == 0) {
            // skip if no neighbors
            continue;
          }

          Real v_out = page_rank[v_idx] / degrees[v_idx];

          auto neighbor_range = *v_it;
          for (auto elt = neighbor_range.begin(); elt != neighbor_range.end(); ++elt) {

            vertex_id_type v = target(G, *elt);

            if (is_same_locality(this_locality_id, G, v)) {
              // handle things locally
              accumulated_contributions[v] = v_out + accumulated_contributions[v];
            }
            else {
              // send to elt's locality
              v_targets.push_back(v);
            }
          }

          // Store any outgoing contributions for this vertex
          if (!v_targets.empty()) {
            for (auto v : v_targets) {
              auto id = vertex_locality(G, v);
              outgoing_packets[id].push_back(std::make_tuple(v, v_out));
            }
          }
        }

        // Now send the outgoing contributions to the other localities
        remote_ops.reserve(remote_ops.size() + outgoing_packets.size());

        page_rank_action_b<Graph, Real> act;
        for (auto&& [id, targets] : outgoing_packets) {
          remote_ops.push_back(hpx::async(act, id, hpx::ref(G), std::move(targets),
                                          hpx::ref(accumulated_contributions)));
        }


        if (!remote_ops.empty()) {
          // wait for all remote operations to finish
          hpx::wait_all(remote_ops);
        }

        return 0;
      }

      template <typename ExPolicy, typename Graph, typename IterB, typename IterE>
      static size_t parallel(ExPolicy&& policy, Graph G, IterB first, IterE last) {
        return 0;
      }
    };

    ////////////////////////////////////////////////////////////////////////////
    template <typename Real>
    struct page_rank_reduce : hpx::parallel::detail::algorithm<page_rank_reduce<Real>, Real> {

      constexpr page_rank_reduce() noexcept
        : hpx::parallel::detail::algorithm<page_rank_reduce, Real>("page_rank_reduce") {}

      template <typename ExPolicy, typename Graph>
      static Real sequential(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index,
                             hpx::partitioned_vector<typename Real> page_rank,
                             hpx::partitioned_vector<typename Real> accumulated_contributions,
                             const Real base_score, const Real damping_factor) {
        // Update page_rank using the accumulated contributions:
        // page_rank[i] += damping_factor * accumulated_contributions[i]

        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;

        std::uint32_t this_locality_id = hpx::get_locality_id();

        using vertex_id_type = typename Graph::vertex_id_type;

        Real local_error = 0.0;

        // for each v in G do
        for (auto v_it = first; v_it != last; ++v_it) {


          size_t v_idx = first_index + std::distance(first, v_it);

          Real old_pr = page_rank[v_idx];
          page_rank[v_idx] = base_score + damping_factor * accumulated_contributions[v_idx];

          // clear the accumulated contributions
          accumulated_contributions[v_idx] = 0.0;

          local_error += fabs(page_rank[v_idx] - old_pr);
        }

        return local_error;
      }

      template <typename ExPolicy, typename Graph>
      static Real parallel(ExPolicy&& policy, Graph G, const size_t first_index,
                           const size_t last_index) {
        return 0;
      }
    };

    /// \endcond
  } // namespace detail


  template <adjacency_list_graph Graph, typename Real = double>
  void partitioned_page_rank_b(Graph& G,
                               hpx::partitioned_vector<typename Graph::vertex_id_type>& degrees,
                               hpx::partitioned_vector<Real>& page_rank,
                               const Real damping_factor = 0.85, const Real threshold = 1.e-4,
                               const size_t max_iters = std::numeric_limits<unsigned int>::max()) {

    auto sizes = G.indices_.get_partition_sizes();

    // hpx::partitioned_vector<vertex_id_type> degrees(
    //   size, std::numeric_limits<vertex_id_type>::max(),
    //   hpx::explicit_container_layout(sizes, G.indices_.get_partition_localities()));
    // degrees.register_as("degrees");

    const Real init_score = 1.0 / page_rank.size();
    const Real base_score = (1.0 - damping_factor) / page_rank.size();

    // Initialize page rank
    hpx::fill(hpx::execution::seq, page_rank.begin(), page_rank.end(), init_score);

    hpx::partitioned_vector<Real> page_rank_accum(
      G.indices_.size(), 0.0,
      hpx::explicit_container_layout(sizes, G.indices_.get_partition_localities()));
    page_rank_accum.register_as("page_rank_accum");

    for (size_t i = 0; i < max_iters; i++) {

      partitioned_algorithm<detail::page_rank_0>(
        hpx::execution::seq, G, hpx::ref(page_rank), hpx::ref(page_rank_accum), hpx::ref(degrees));

      auto errors = partitioned_algorithm<detail::page_rank_reduce<Real>>(
        hpx::execution::seq, G, hpx::ref(page_rank), hpx::ref(page_rank_accum), base_score,
        damping_factor);

      auto error = std::transform_reduce(
        errors.begin(), errors.end(), Real(0.0), [](Real count, Real curr) { return count + curr; },
        [](auto&& f) { return f.get(); });

      if (error < threshold)
        break;
    }
  }
} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_PAGE_RANK_B_HPP
