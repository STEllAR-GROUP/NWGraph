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

#ifndef NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_1_HPP
#define NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_1_HPP

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

    //    template <typename Graph, typename Vertex>
    //    hpx::id_type vertex_locality(Graph const& G, Vertex v) {
    //      using traits = hpx::traits::segmented_iterator_traits<decltype(G.begin().index())>;
    //      return traits::get_id(traits::segment((G.begin() + v).index()));
    //    }
    //
    //    template <typename Graph, typename Vertex>
    //    bool is_same_locality(std::uint32_t this_locality_id, Graph const& G, Vertex v) {
    //      return this_locality_id == hpx::naming::get_locality_id_from_id(vertex_locality(G, v));
    //    }

    ////////////////////////////////////////////////////////////////////////////
    // handle counting of triangles on target locality
    template <typename Graph>
    static size_t triangle_counter_1(
      Graph G,
      std::vector<std::tuple<typename Graph::vertex_id_type,
                             std::vector<std::tuple<typename Graph::vertex_id_type>>>> const&
        targets) {

      size_t triangles = 0;
      for (auto&& [v, neighbors] : targets) {
        triangles += nw::graph::intersection_size(neighbors, G[v]);
      }
      return triangles;
    }

    template <typename Graph>
    struct triangle_count_action_1
      : hpx::actions::action<decltype(&triangle_counter_1<Graph>), &triangle_counter_1<Graph>,
                             triangle_count_action_1<Graph>> {};

    ////////////////////////////////////////////////////////////////////////////
    struct triangle_count_1 : hpx::parallel::detail::algorithm<triangle_count_1, size_t> {

      // triangle counting driver for one of the partitions
      constexpr triangle_count_1() noexcept
        : hpx::parallel::detail::algorithm<triangle_count_1, size_t>("triangle_count_1") {}

      template <typename ExPolicy, typename Graph>
      static size_t sequential(ExPolicy&& policy, Graph G, size_t first_index, size_t last_index) {

        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;

        size_t triangles = 0;
        std::uint32_t this_locality_id = hpx::get_locality_id();

        using vertex_id_type = typename Graph::vertex_id_type;
        std::vector<std::tuple<vertex_id_type>> neighbors;
        std::map<hpx::id_type,
                 std::vector<std::tuple<vertex_id_type, std::vector<std::tuple<vertex_id_type>>>>>
          remote_counts;
        std::vector<vertex_id_type> v_targets;

        // for each v in G do
        for (auto v_it = first; v_it != last; ++v_it) {

          v_targets.resize(0);
          neighbors.resize(0);

          auto neighbor_range = *v_it;
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
            for (auto v : v_targets) {
              auto id = vertex_locality(G, v);
              remote_counts[id].push_back(std::make_tuple(v, neighbors));
            }
          }
        }

        std::vector<hpx::future<size_t>> counts;

        triangle_count_action_1<Graph> act;
        for (auto&& [id, targets] : remote_counts) {
          counts.push_back(hpx::async(act, hpx::colocated(id), hpx::ref(G), std::move(targets)));
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
  size_t partitioned_triangle_count_1(Graph& G) {
    auto counts = partitioned_algorithm<detail::triangle_count_1>(hpx::execution::seq, G);
    return std::transform_reduce(
      counts.begin(), counts.end(), size_t(0),
      [](size_t count, size_t curr) { return count + curr; }, [](auto&& f) { return f.get(); });
  }
} // namespace nw::graph

#endif //  NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_1_HPP
