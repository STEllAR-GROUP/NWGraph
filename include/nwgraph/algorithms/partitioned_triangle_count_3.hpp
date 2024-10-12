/**
 * @file triangle_count.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * @authors
 *   Andrew Lumsdaine
 *   Tony Liu
 *   Kevin Deweese
 *
 */

#ifndef NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_3_HPP
#define NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_3_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/algorithms/partitioned_algorithm.hpp"
#include "nwgraph/algorithms/triangle_count.hpp"

#include <algorithm>
#include <cstddef>
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
    struct triangle_count_3 : hpx::parallel::detail::algorithm<triangle_count_3, size_t> {

      // triangle counting driver for one of the partitions
      constexpr triangle_count_3() noexcept
        : hpx::parallel::detail::algorithm<triangle_count_3, size_t>("triangle_count_3") {}

      template <typename ExPolicy, typename Graph>
      static size_t sequential(ExPolicy&&, Graph G, size_t first_index, size_t last_index,
                               size_t batchsize) {

        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;

        std::uint32_t this_locality_id = hpx::get_locality_id();

        using vertex_id_type = typename Graph::vertex_id_type;
        using target_list_t =
          std::vector<std::tuple<vertex_id_type, std::vector<std::tuple<vertex_id_type>>>>;
        using remote_counts_t = std::map<hpx::id_type, target_list_t>;

        // use half of the available cores for parallelizing the loop
        auto p = hpx::execution::par;
        size_t const cores = hpx::parallel::execution::processing_units_count(
          p.parameters(), p.executor(), hpx::chrono::null_duration, last_index - first_index);
        auto policy = hpx::parallel::execution::with_processing_units_count(
          p, (std::max)(cores / 2, size_t(1)));

        // for each v in G do
        safe_object<std::tuple<size_t, remote_counts_t, std::vector<hpx::future<size_t>>>>
          remote_counts;

        auto tc = [&](auto&& neighbor_range)
        {
          size_t triangles = 0;
          std::vector<vertex_id_type> v_targets;
          std::vector<std::tuple<vertex_id_type>> neighbors;

          for (auto elt = neighbor_range.begin(); elt != neighbor_range.end(); ++elt) {

            vertex_id_type v = target(G, *elt);

            if (is_same_locality(this_locality_id, G, v)) {
              // handle things locally
              triangles += nw::graph::intersection_size(neighbor_range, G[v]);
            }
            else {
              // send our neighbors to elt's locality
              v_targets.push_back(v);
            }

            // collect all neighbor vertex ids for v_it
            neighbors.push_back(std::make_tuple(v));
          }

          // launch the remote operations for the current vertex (if any)
          if (!v_targets.empty()) {
            auto& rc = remote_counts.get();

            for (auto v : v_targets) {
              auto id = vertex_locality(G, v);
              auto& targets = std::get<1>(rc)[id];

              // if batch size has been reached, trigger async operation
              if (targets.size() >= batchsize) {
                triangle_count_action<Graph> act;
                std::get<2>(rc).push_back(hpx::async(act, id, hpx::ref(G), std::move(targets)));
                targets = target_list_t{};
              }

              targets.push_back(std::make_tuple(v, neighbors));
            }
          }
          std::get<0>(remote_counts.get()) += triangles;
        };
        hpx::for_each(policy, first, last, tc);

        // collect remote counts
        std::vector<hpx::future<size_t>> counts;
        size_t triangles = 0;

        remote_counts.reduce(
          [&](auto&& data)
          {
            // accumulate all local triangle counts
            triangles += std::get<0>(data);

            // send remaining pending messages
            for (auto&& [id, targets] : std::get<1>(data)) {
              triangle_count_action<Graph> act;
              counts.push_back(hpx::async(act, id, hpx::ref(G), std::move(targets)));
            }

            std::move(std::get<2>(data).begin(), std::get<2>(data).end(),
                      std::back_inserter(counts));
          });

        if (!counts.empty()) {
          // wait for all remote operations to finish
          hpx::wait_all(counts);
          return std::transform_reduce(
            counts.begin(), counts.end(), triangles, [](size_t count, size_t curr)
            { return count + curr; }, [](auto&& f) { return f.get(); });
        }

        return triangles;
      }

      template <typename ExPolicy, typename Graph>
      static size_t parallel(ExPolicy&& policy, Graph const& G, size_t first_index,
                             size_t last_index, size_t batchsize) {
        return 0;
      }
    };
    /// \endcond
  } // namespace detail

  /**
   * @brief Two-dimensional triangle counting, parallel version.
   *
   * @tparam Graph adjacency_list_graph
   * @param G graph
   * @return size_t number of triangles
   */

  template <adjacency_list_graph Graph>
  size_t partitioned_triangle_count_3(Graph& G, size_t batchsize) {
    auto counts =
      partitioned_algorithm<detail::triangle_count_3>(hpx::execution::seq, G, batchsize);
    return std::transform_reduce(
      counts.begin(), counts.end(), size_t(0),
      [](size_t count, size_t curr) { return count + curr; }, [](auto&& f) { return f.get(); });
  }
} // namespace nw::graph

#endif //  NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_3_HPP
