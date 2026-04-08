/**
 * @file triangle_count.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * @authors
 *   Hartmut Kaiser
 *
 */

#ifndef NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_0_HPP
#define NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_0_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/distributed/algorithms/algorithm.hpp"
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

    namespace triangle_count_0_impl {

    ////////////////////////////////////////////////////////////////////////////
    // handle counting of triangles on target locality
    template <typename Graph>
    static size_t
    triangle_counter_0(Graph G, typename Graph::vertex_id_type v,
                       std::vector<std::tuple<typename Graph::vertex_id_type>> const& neighbors) {
      using vertex_id_type = typename Graph::vertex_id_type;
      auto cmp = [&G](auto&& lhs, auto&& rhs)
      { return static_cast<vertex_id_type>(target(G, lhs)) < static_cast<vertex_id_type>(target(G, rhs)); };
      auto target_neighbors = G[v];
      return nw::graph::intersection_size(neighbors, target_neighbors, cmp);
    }

    template <typename Graph>
    struct triangle_count_action_0
      : hpx::actions::action<decltype(&triangle_counter_0<Graph>), &triangle_counter_0<Graph>,
                             triangle_count_action_0<Graph>> {};

    ////////////////////////////////////////////////////////////////////////////
    struct triangle_count_0 : hpx::parallel::detail::algorithm<triangle_count_0, size_t> {

      // triangle counting driver for one of the partitions
      constexpr triangle_count_0() noexcept
        : hpx::parallel::detail::algorithm<triangle_count_0, size_t>("triangle_count_0") {}

      template <typename ExPolicy, typename Graph>
      static size_t sequential(ExPolicy&& policy, Graph G, size_t first_index, size_t last_index) {

        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;

        size_t triangles = 0;
        std::uint32_t this_locality_id = hpx::get_locality_id();

        using vertex_id_type = typename Graph::vertex_id_type;
        std::vector<hpx::future<size_t>> counts;
        std::vector<std::tuple<vertex_id_type>> neighbors;
        std::vector<vertex_id_type> targets;
        auto cmp = [&G](auto&& lhs, auto&& rhs)
        { return static_cast<vertex_id_type>(target(G, lhs)) < static_cast<vertex_id_type>(target(G, rhs)); };

        // for each v in G do
        for (auto v_it = first; v_it != last; ++v_it) {

          targets.resize(0);
          neighbors.resize(0);

          auto neighbor_range = *v_it;
          for (auto elt = neighbor_range.begin(); elt != neighbor_range.end(); ++elt) {

            vertex_id_type v = static_cast<vertex_id_type>(target(G, *elt));

            if (nw::graph::detail::is_same_locality(this_locality_id, G, v)) {
              // handle things locally
              auto target_neighbors = G[v];
              triangles += nw::graph::intersection_size(neighbor_range, target_neighbors, cmp);
            }
            else {
              // send our neighbors to elt's locality
              targets.push_back(v);
            }

            // collect all neighbor vertex ids for v_it
            neighbors.push_back(std::make_tuple(v));
          }

          // launch the remote operations for the current vertex (if any)
          if (!targets.empty()) {
            counts.reserve(counts.size() + targets.size());

            triangle_count_action_0<Graph> act;
            for (auto v : targets) {
              auto id = nw::graph::detail::vertex_locality(G, v);
              counts.push_back(hpx::async(act, id, hpx::ref(G), v, neighbors));
            }
          }
        }

        if (!counts.empty()) {
          // wait for all remote operations to finish
          hpx::wait_all(counts);
          return std::transform_reduce(
            counts.begin(), counts.end(), triangles, [](size_t count, size_t curr)
            { return count + curr; }, [](auto&& f) { return f.get(); });
        }

        return triangles;
      }

      template <typename ExPolicy, typename Graph, typename IterB, typename IterE>
      static size_t parallel(ExPolicy&& policy, Graph const& G, IterB first, IterE last) {
        return 0;
      }
    };

    } // namespace triangle_count_0_impl
    /// \endcond
  } // namespace detail

  /**
   * @brief Two-dimensional triangle counting, parallel version.
   *
   * @tparam Graph adjacency_list_graph
   * @param G graph
   * @param threads number of threads
   * @return size_t number of triangles
   */
  template <adjacency_list_graph Graph>
  size_t partitioned_triangle_count_0(Graph& G) {
    auto counts = partitioned_algorithm<detail::triangle_count_0_impl::triangle_count_0>(
      hpx::execution::seq, G);
    return std::transform_reduce(
      counts.begin(), counts.end(), size_t(0),
      [](size_t count, size_t curr) { return count + curr; }, [](auto&& f) { return f.get(); });
  }
} // namespace nw::graph

#endif //  NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_0_HPP
