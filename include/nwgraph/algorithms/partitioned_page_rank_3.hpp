/**
 * @file partitioned_page_rank_3.hpp
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

#ifndef NW_GRAPH_PARTITIONED_page_rank_3_HPP
#define NW_GRAPH_PARTITIONED_page_rank_3_HPP

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

///////////////////////////////////////////////////////////////////////////////
namespace nw::graph {

  namespace detail {


    template <typename Graph, typename Real>
    using remote_results_type = std::vector<std::tuple<typename Graph::vertex_id_type, Real>>;

    template <typename Graph, typename Real>
    static remote_results_type<Graph, Real> do_page_rank_packet_3(
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
    struct page_rank_action_3
      : hpx::actions::action<decltype(&do_page_rank_packet_3<Graph, Real>),
                             &do_page_rank_packet_3<Graph, Real>, page_rank_action_3<Graph, Real>> {
    };

    ////////////////////////////////////////////////////////////////////////////
    template <typename Real>
    struct page_rank_3 : hpx::parallel::detail::algorithm<page_rank_3<Real>, Real> {

      // page rank driver for one of the partitions
      constexpr page_rank_3() noexcept
        : hpx::parallel::detail::algorithm<page_rank_3, Real>("page_rank_3") {}


      template <typename ExPolicy, typename Graph>
      static Real sequential(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index, hpx::partitioned_vector<Real> page_rank,
                             hpx::partitioned_vector<Real> accumulated_contributions,
                             hpx::partitioned_vector<typename Graph::vertex_id_type> degrees,
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

        hpx::id_type this_locality = hpx::find_here();


        using vertex_id_type = typename Graph::vertex_id_type;

        using idx_t = typename Graph::vertex_id_type;

        // A packet is a tuple of (destination vertex, query vertex)
        using packet_vec_t = std::vector<std::tuple<idx_t, idx_t>>;

        // Seperate outgoing packet vectors as to which locality they need to be sent
        using locality_map_t = std::map<hpx::id_type, packet_vec_t>;

        // Let each thread work on its own
        safe_object<std::tuple<locality_map_t, std::vector<hpx::future<void>>>> thread_locals;

        // This function works on a single partition
        // Lets query only once for the local iterators
        auto ac_first = accumulated_contributions.get_local_iterator(first_index).local();
        auto pr_first = page_rank.get_local_iterator(first_index).local();
        auto deg_first = degrees.get_local_iterator(first_index).local();


        auto send_packets = [&](hpx::id_type id, packet_vec_t&& packet)
        {
          using action_t = page_rank_action_3<Graph, Real>;
          auto f =
            hpx::async(action_t(), id, std::move(packet), hpx::ref(page_rank), hpx::ref(degrees));
          // When the result has been received, add the contributions
          return f.then(
            [&](auto&& f)
            {
              auto incoming_packet = f.get();
              for (auto&& [v, incoming_val] : incoming_packet) {
                auto offs = v - first_index;
                auto ac_iter = ac_first + offs;
                std::atomic_ref(*ac_iter) += incoming_val;
              }
            });
        };


        // Thread local copies of the iterators
        auto tl_iters = safe_object(std::make_tuple(ac_first, pr_first, deg_first));

        // Iterate over all local vertices
        hpx::experimental::for_loop(hpx::execution::par, first_index, last_index,
            [&](auto&& u) {
   
            //auto ac_iter = accumulated_contributions.get_local_iterator(u).local();

            // For each neighbor of u
            for (auto&& elt : G[u]) {
                
                idx_t v = target(G, elt);
                auto v_locality = vertex_locality(G, v);
                if (v_locality == this_locality) {
                  auto& ac_iter = std::get<0>(tl_iters.get()); //+ (u - first_index);
                  ac_iter += (u - first_index);
                  auto offs = v - first_index;
                  auto& pr_iter = std::get<1>(tl_iters.get()); //+ offs;
                  auto& deg_iter = std::get<2>(tl_iters.get()); //+ offs;
                  pr_iter += offs;
                  deg_iter += offs;
                  //auto pr_iter = page_rank.get_local_iterator(v).local();
                  //auto deg_iter = degrees.get_local_iterator(v).local();
                  *(ac_iter) += *(pr_iter) / *(deg_iter);
                  ac_iter -= (u - first_index);
                  pr_iter -= offs;
                  deg_iter -= offs;
                } else{
                  // Attributes of v are remote, so page_rank[v] / degrees[v] needs to be computed
                  // remotely and sent back to the current 
                  auto& [map, futures] = thread_locals.get();
                  auto& packets = map[v_locality];
                  packets.push_back(std::make_tuple(u, v));
                  // If packet for that locality is big enough, send it right away
                  if (packets.size() >= batchsize) {
                    auto tmp = std::move(packets);
                    packets = packet_vec_t{};
                    auto fut = send_packets(v_locality, std::move(tmp));
                    futures.push_back(std::move(fut));
                   }
                }
            }
            });

        // Flush remaining packets and wait for remote results to be received
        std::vector<hpx::future<void>> all_futures;

        thread_locals.reduce(
            [&](auto&& data)
          {
            auto& map = std::get<0>(data);
            for (auto&& [id, packets] : map) {
              if (!packets.empty()) {
                auto tmp = std::move(packets);
                packets = packet_vec_t{};
                all_futures.push_back(send_packets(id, std::move(tmp)));
              }
            }

            auto& futures = std::get<1>(data);  
            std::move(futures.begin(), futures.end(), std::back_inserter(all_futures));
          });

        hpx::wait_all(all_futures);

        // (local) accumulated result is now fully computed
        // update the local page rank
        Real local_error = 0.0;

        hpx::experimental::for_loop(hpx::execution::par, first_index, last_index, 
            hpx::experimental::reduction_plus(local_error),
            [&](auto&& u, auto& error)
            {
                auto offs = u - first_index;
                auto& acc_iter = std::get<0>(tl_iters.get()); // +offs;
                acc_iter += offs;
                //auto acc_iter = accumulated_contributions.get_local_iterator(u).local();
                Real z = *acc_iter;

                auto& pr_iter = std::get<1>(tl_iters.get());// + offs;
                pr_iter += offs;
                //auto pr_iter = page_rank.get_local_iterator(u).local();
                auto old_rank = *pr_iter;
                auto new_rank = base_score + damping_factor * z;
                *pr_iter = new_rank;
                error += fabs(new_rank - old_rank);
                *acc_iter = 0.0;

                acc_iter -= offs;
                pr_iter -= offs;
            });

        return local_error;
      }

      template <typename ExPolicy, typename Graph, typename IterB, typename IterE>
      static size_t parallel(ExPolicy&& policy, Graph G, IterB first, IterE last, size_t batchsize) {
        return 0;
      }
    };


    /// \endcond
  } // namespace detail

  template <adjacency_list_graph Graph, typename Real = double>
  void partitioned_page_rank_3(Graph& G,
                               hpx::partitioned_vector<typename Graph::vertex_id_type>& degrees,
                               hpx::partitioned_vector<Real>& page_rank,
                               const Real damping_factor = 0.85, const Real threshold = 1.e-4,
                               const size_t max_iters = std::numeric_limits<unsigned int>::max(),
                               size_t batchsize = 1000) {

    const Real init_score = 1.0 / G.size();
    const Real base_score = (1.0 - damping_factor) / G.size();

    // Initialize page ranks
    hpx::fill(hpx::execution::seq, page_rank.begin(), page_rank.end(), init_score);


    // Create buffer vector needed for algorithm
    hpx::partitioned_vector<Real> accum(
      page_rank.size(),
      hpx::explicit_container_layout(page_rank.get_partition_sizes(),
                                     page_rank.get_partition_localities()));
    accum.register_as("accum");

    // Initialize accumulations
    hpx::fill(hpx::execution::seq, accum.begin(), accum.end(), 0.0f);


    // Run PR iterations
    std::vector<hpx::future<Real>> errors;


    for (size_t iter = 0; iter < max_iters; ++iter) {
      std::cout << "----- Iteration " << iter << " ----- ";

      // Keep time
      auto t_start = std::chrono::high_resolution_clock::now();

      //for (auto i = 0; i < page_rank.size(); ++i) {
      //  std::cout << "Node " << i << " : " << page_rank[i] << std::endl;
      //}

      errors = partitioned_segmented_algorithm<detail::page_rank_3<Real>>(
        hpx::execution::seq, G, hpx::ref(page_rank), hpx::ref(accum), hpx::ref(degrees), base_score,
        damping_factor, batchsize);

      Real error = std::transform_reduce(
        errors.begin(), errors.end(), Real(0.0), [](Real count, Real curr) { return count + curr; },
        [](auto&& f) { return f.get(); });

      auto t_end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = t_end - t_start;
      std::cout << " (dt = " << diff.count() << " sec) " << std::endl;

      if (error < threshold)
        break;
    }
  }
} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_page_rank_3_HPP
