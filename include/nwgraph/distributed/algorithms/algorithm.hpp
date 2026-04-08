/**
 * @file partitioned_algorithm.hpp
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

#ifndef NW_GRAPH_PARTITIONED_ALGORITHM_HPP
#define NW_GRAPH_PARTITIONED_ALGORITHM_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/algorithms/triangle_count.hpp"
#include "nwgraph/distributed/graph.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <ranges>
#include <utility>
#include <type_traits>
#include <vector>

#include <hpx/async_combinators/wait_all.hpp>
#include <hpx/executors/execution_policy.hpp>
#include <hpx/include/partitioned_vector_predef.hpp>
#include <hpx/parallel/segmented_algorithms/detail/dispatch.hpp>
#include <hpx/parallel/util/detail/algorithm_result.hpp>
#include <hpx/concurrency/cache_line_data.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace nw::graph {

  namespace detail {
    template <typename Graph>
    decltype(auto) remote_graph_ref(Graph& G) {
      return hpx::ref(G);
    }

    template <typename Graph, typename Vertex>
    hpx::id_type vertex_locality(Graph const& G, Vertex v) {
      using graph_type = std::remove_reference_t<Graph>;
      using vertex_id_type = vertex_id_t<graph_type>;
      return static_cast<hpx::id_type>(vertex_partition(G, static_cast<vertex_id_type>(v)));
    }

    template <typename Graph, typename Vertex>
    std::uint32_t vertex_locality_id(Graph const& G, Vertex v) {
      return hpx::naming::get_locality_id_from_id(vertex_locality(G, v));
    }

    template <typename Graph, typename Vertex>
    bool is_same_locality(std::uint32_t this_locality_id, Graph const& G, Vertex v) {
      return this_locality_id == hpx::naming::get_locality_id_from_id(vertex_locality(G, v));
    }

    ////////////////////////////////////////////////////////////////////////////
    template <typename T>
    struct safe_object {
      safe_object()
        : data_(hpx::get_os_thread_count()) {}

      safe_object(T const& init)
        : data_(hpx::get_os_thread_count(), init) {}

      safe_object(safe_object const& rhs) = delete;
      safe_object(safe_object&& rhs) noexcept = default;

      safe_object& operator=(safe_object const& rhs) = delete;
      safe_object& operator=(safe_object&& rhs) noexcept = default;

      T& get() { return data_[hpx::get_worker_thread_num()].data_; }

      T const& get() const { return data_[hpx::get_worker_thread_num()].data_; }

      template <typename F>
      void reduce(F const& f) {
        for (auto&& d : std::move(data_)) {
          f(std::move(d.data_));
        }
      }

    private:
      std::vector<hpx::util::cache_line_data<T>> data_;
    };

    ////////////////////////////////////////////////////////////////////////////
    // handle counting of triangles on target locality
    template <typename Graph>
    static size_t triangle_counter(
      Graph G,
      std::vector<std::tuple<vertex_id_t<std::remove_reference_t<Graph>>,
                             std::vector<std::tuple<vertex_id_t<std::remove_reference_t<Graph>>>>>> const&
        targets) {

      size_t triangles = 0;
      for (auto&& [v, neighbors] : targets) {
        triangles += nw::graph::intersection_size(neighbors, G[v]);
      }
      return triangles;
    }

    template <typename Graph>
    struct triangle_count_action
      : hpx::actions::action<decltype(&triangle_counter<Graph>), &triangle_counter<Graph>,
                             triangle_count_action<Graph>> {};

    ////////////////////////////////////////////////////////////////////////////
    // handle counting of triangles on target locality
    template <typename Graph>
    static size_t triangle_counter_1(
      Graph G,
      std::vector<std::tuple<std::vector<vertex_id_t<std::remove_reference_t<Graph>>>,
                             std::vector<std::tuple<vertex_id_t<std::remove_reference_t<Graph>>>>>> const&
        targets) {

      size_t triangles = 0;
      for (auto&& target : targets) {
        for (auto v : std::get<0>(target)) {
          triangles += nw::graph::intersection_size(std::get<1>(target), G[v]);
        }
      }
      return triangles;
    }

    template <typename Graph>
    struct triangle_count_action_1
      : hpx::actions::action<decltype(&triangle_counter_1<Graph>), &triangle_counter_1<Graph>,
                             triangle_count_action_1<Graph>> {};

  } // namespace detail

  template <typename G>
  concept partitioned_algorithm_graph = partitioned_graph<G>;

  /**
   * @brief Generic segmented distributed algorithm implementation.
   *
   * This entry point dispatches once per concrete graph segment and passes the
   * segment's global first/last indices to the algorithm implementation.
   */
  template <typename Algorithm, typename ExPolicy, partitioned_algorithm_graph Graph,
            typename... Ts>
  [[gnu::noinline]] auto partitioned_algorithm(ExPolicy&& policy, Graph& G, Ts&&... ts) {
    using algorithm_t = Algorithm;
    using result_t = typename algorithm_t::result_type;

    std::vector<hpx::future<result_t>> results;
    auto segments = partition_segments(G);
    results.reserve(segments.size());

    for (auto const& segment : segments) {
      results.push_back(dispatch_async(
        segment.locality, algorithm_t{}, policy, std::true_type(),
        detail::remote_graph_ref(G),
        segment.first, segment.last, ts...));
    }

    hpx::wait_all(results);
    return results;
  }
} // namespace nw::graph

#endif //  NW_GRAPH_PARTITIONED_ALGORITHM_HPP
