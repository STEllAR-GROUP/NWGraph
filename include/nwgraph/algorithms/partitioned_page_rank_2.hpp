/**
 * @file partitioned_page_rank_2.hpp
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

#ifndef NW_GRAPH_PARTITIONED_PAGE_RANK_2_HPP
#define NW_GRAPH_PARTITIONED_PAGE_RANK_2_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/algorithms/partitioned_algorithm.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <vector>

#include <hpx/async_combinators/wait_all.hpp>
#include <hpx/executors/execution_policy.hpp>
#include <hpx/include/partitioned_vector_predef.hpp>
#include <hpx/parallel/segmented_algorithms/detail/dispatch.hpp>
#include <hpx/parallel/util/detail/algorithm_result.hpp>

#include <hpx/threading_base/annotated_function.hpp>
#include <hpx/threading_base/scoped_annotation.hpp>

#include <hpx/barrier.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace nw::graph {

  namespace detail {

    template <typename Graph, typename Real>
    static void do_page_rank_packet_2(
      std::vector<std::tuple<typename Graph::vertex_id_type, Real>>&& incoming_packet,
      hpx::partitioned_vector<Real> page_rank) {

      hpx::scoped_annotation annotation("do_page_rank_packet_2");

      for (auto&& [v_dest, val] : incoming_packet) {
        auto pr_iter = page_rank.get_local_iterator(v_dest).local();
        *pr_iter += val;
      }
    }

    template <typename Graph, typename Real>
    struct page_rank_action_2
      : hpx::actions::action<decltype(&do_page_rank_packet_2<Graph, Real>),
                             &do_page_rank_packet_2<Graph, Real>, page_rank_action_2<Graph, Real>> {
    };

    ////////////////////////////////////////////////////////////////////////////
    template <typename Real>
    struct page_rank_2 : hpx::parallel::detail::algorithm<page_rank_2<Real>, Real> {

      // page rank driver for one of the partitions
      constexpr page_rank_2() noexcept
        : hpx::parallel::detail::algorithm<page_rank_2, Real>("page_rank_2") {}


      template <typename ExPolicy, typename Graph>
      static Real sequential(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index, hpx::partitioned_vector<Real> from_page_rank,
                             hpx::partitioned_vector<Real> to_page_rank, Real base_score,
                             Real damping_factor) {

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


        //


        hpx::id_type this_locality_id = hpx::find_here();


        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;


        using vertex_id_type = typename Graph::vertex_id_type;

        auto to_pr_iter = to_page_rank.get_local_iterator(first_index).local();
        {
          hpx::scoped_annotation annotation("page_rank_2::init_base_score");
          // Initialize values to base_score
          for (auto v_it = first; v_it != last; v_it++, to_pr_iter++) {
            *to_pr_iter = base_score;
          }
        }

        //
        std::vector<vertex_id_type> v_remote;

        std::map<hpx::id_type, std::vector<std::tuple<vertex_id_type, Real>>> outgoing_packets;
        std::vector<hpx::future<void>> remote_ops;
        auto from_pr_iter = from_page_rank.get_local_iterator(first_index).local();

        {
          hpx::scoped_annotation annotation("page_rank_2::compute_contributions");
          for (auto v_it = first; v_it != last; v_it++, from_pr_iter++) {
            auto v_id = v_it.index();
            auto neighbour_range = *v_it;
            auto out_degree = neighbour_range.size();
            auto out_rank = damping_factor * (*from_pr_iter) / out_degree;
            // for each neighbour
            for (auto&& e : neighbour_range) {
              auto v = target(G, e);
              auto loc_id = vertex_locality(G, v);
              if (loc_id == this_locality_id) {
                // local contribution
                auto local_it = to_page_rank.get_local_iterator(v).local();
                *local_it += out_rank;
              }
              else {
                // remote contribution
                outgoing_packets[loc_id].push_back({v, out_rank});
              }
            }
          }
        }

        {
          hpx::scoped_annotation annotation("page_rank_2::remote_contributions");
          // Send remote contributions
          for (auto&& [id, packet] : outgoing_packets) {
            using act_t = page_rank_action_2<Graph, Real>;
            remote_ops.push_back(hpx::async<act_t>(id, std::move(packet), hpx::ref(to_page_rank)));
          }

          // Wait for remote contributions to be received
          hpx::wait_all(remote_ops);
        }


        // Barrier, only return after all localities have finished their work
        // TODO: limit to localities that own part of the distributed vector
        //
        // Hack, assume this function runs once for each graph partition
        auto get_segment_barrier = [](auto& partitioned_vec, size_t idx)
        {
          auto n_parts = partitioned_vec.partitions().size();
          auto part_0 = partitioned_vec.get_segment_iterator(0);
          auto part_iter = partitioned_vec.get_segment_iterator(idx);
          auto part_idx = part_iter - part_0;
          return hpx::distributed::barrier("nw::graph::partitioned_page_rank_2_barrier", n_parts,
                                           part_idx);
        };

        {
          hpx::scoped_annotation annotation("page_rank_2::barrier");
          auto barrier = get_segment_barrier(to_page_rank, first_index);
          barrier.wait();
        }


        // (local) accumulated result is now fully computed
        // Calculate the error
        Real local_error = 0.0;

        to_pr_iter = to_page_rank.get_local_iterator(first_index).local();
        from_pr_iter = from_page_rank.get_local_iterator(first_index).local();

        {
          hpx::scoped_annotation annotation("page_rank_2::compute_error");
          for (auto v_it = first; v_it != last; ++v_it, ++to_pr_iter, ++from_pr_iter) {
            local_error += fabs(*to_pr_iter - *from_pr_iter);
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
  void partitioned_page_rank_2(Graph& G,
                               hpx::partitioned_vector<Real>& page_rank,
                               const Real damping_factor = 0.85, const Real threshold = 1.e-4,
                               const size_t max_iters = std::numeric_limits<unsigned int>::max()) {

    hpx::scoped_annotation annotation("partitioned_page_rank_2");

    const Real init_score = 1.0 / G.size();
    const Real base_score = (1.0 - damping_factor) / G.size();

    // Initialize page ranks
    hpx::fill(hpx::execution::seq, page_rank.begin(), page_rank.end(), init_score);


    // Create buffer vector needed for algorithm
    hpx::partitioned_vector<Real> page_rank_2(
      page_rank.size(),
      hpx::explicit_container_layout(page_rank.get_partition_sizes(),
                                     page_rank.get_partition_localities()));
    page_rank_2.register_as("page_rank_2");


    // Run PR iterations
    std::vector<hpx::future<Real>> errors;
    bool to_page_rank_2 = false;

    auto page_rank_iter = [&](auto&& from_pr, auto&& to_pr)
    {
      return partitioned_algorithm<detail::page_rank_2<Real>>(
        hpx::execution::seq, G, hpx::ref(from_pr), hpx::ref(to_pr), base_score, damping_factor);
    };

    for (size_t iter = 0; iter < max_iters; ++iter) {
      std::cout << "----- Iteration " << iter << " ----- " << std::endl;
      hpx::scoped_annotation annotation("page_rank_iter");

      //for (auto i = 0; i < page_rank.size(); ++i) {
      //  std::cout << "Node " << i << " : " << page_rank[i] << std::endl;
      //}

      if (to_page_rank_2) {
        errors = page_rank_iter(hpx::ref(page_rank), hpx::ref(page_rank_2));
        to_page_rank_2 = false;
      }
      else {
        errors = page_rank_iter(hpx::ref(page_rank_2), hpx::ref(page_rank));
        to_page_rank_2 = true;
      }

      Real error = std::transform_reduce(
        errors.begin(), errors.end(), Real(0.0), [](Real count, Real curr) { return count + curr; },
        [](auto&& f) { return f.get(); });

      if (error < threshold)
        break;
    }

    // Copy latest page rank to result, if needed
    if (!to_page_rank_2) {
      // TODO
    }
  }
} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_PAGE_RANK_0_HPP
