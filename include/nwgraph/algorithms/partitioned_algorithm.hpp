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

#include "nwgraph/partition.hpp"
#include "nwgraph/algorithms/triangle_count.hpp"

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <utility>
#include <type_traits>
#include <vector>

#include <hpx/async_combinators/wait_all.hpp>
#include <hpx/algorithms/traits/is_value_proxy.hpp>
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
      if constexpr (tag_invocable<remote_ref_tag, Graph&>) {
        return remote_ref(G);
      }
      else {
        return hpx::ref(G);
      }
    }

    template <typename Container, typename Partition>
      requires partition_token<std::remove_cvref_t<Partition>>
    decltype(auto) partition_data(Container& container, Partition const& partition) {
      if constexpr (tag_invocable<local_view_tag, Container&, Partition>) {
        return local_view(container, partition);
      }
      else {
        return (container);
      }
    }

    template <typename Row>
    auto iterable_row(Row&& row) {
      using row_type = std::remove_cvref_t<Row>;

      if constexpr (hpx::traits::is_value_proxy_v<row_type>) {
        return static_cast<hpx::traits::proxy_value_t<row_type>>(row);
      }
      else {
        return std::forward<Row>(row);
      }
    }

    template <typename Graph, typename Edge>
    auto edge_target(Graph const& graph, Edge&& edge) {
      if constexpr (requires { target(graph, std::forward<Edge>(edge)); }) {
        return target(graph, std::forward<Edge>(edge));
      }
      else {
        using graph_type = std::remove_reference_t<Graph>;
        return static_cast<vertex_id_t<graph_type>>(std::forward<Edge>(edge));
      }
    }

    template <typename Row>
    auto row_size(Row&& row) {
      auto normalized_row = iterable_row(std::forward<Row>(row));
      return std::ranges::size(normalized_row);
    }

    template <typename Algorithm, typename = void>
    struct is_partition_aware_algorithm : std::false_type {};

    template <typename Algorithm>
    struct is_partition_aware_algorithm<Algorithm, std::void_t<decltype(Algorithm::partition_aware)>>
      : std::bool_constant<Algorithm::partition_aware> {};

    template <typename Graph, typename Vertex>
    hpx::id_type vertex_locality(Graph const& G, Vertex v) {
      return partition_locality(vertex_partition(G, v));
    }

    template <typename Graph, typename Vertex>
    std::uint32_t vertex_locality_id(Graph const& G, Vertex v) {
      return hpx::naming::get_locality_id_from_id(vertex_locality(G, v));
    }

    template <typename Graph, typename Vertex>
    bool is_same_locality(std::uint32_t this_locality_id, Graph const& G, Vertex v) {
      return this_locality_id == hpx::naming::get_locality_id_from_id(vertex_locality(G, v));
    }

    template <typename Traits, typename SegIter, typename LocalIter>
    auto global_index(SegIter const& seg, LocalIter const& local) {

      auto composed = Traits::compose(seg, local);
      return composed.get_global_index();
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

    template <typename G>
    concept legacy_segmented_partitioned_graph =
      adjacency_list_graph<G> &&
      requires(G& g) {
        g.begin().index();
        g.end().index();
      };
  } // namespace detail

  template <typename G>
  concept partitioned_algorithm_graph = partitioned_graph<G>;

  /**
   * @brief Generic partition-aware algorithm implementation.
   *
  * This is the semantic distributed graph entry point. Algorithms dispatched
  * here operate on partition descriptors and global index ranges. `local_view`
  * may still be used as an optional optimization, but is not part of the
  * required graph integration contract.
   */
  template <typename Algorithm, typename ExPolicy, partitioned_graph Graph, typename... Ts>
    requires detail::is_partition_aware_algorithm<Algorithm>::value
  [[gnu::noinline]] auto partitioned_algorithm(ExPolicy&& policy, Graph& G, Ts&&... ts) {
    using algorithm_t = Algorithm;
    using result_t = typename algorithm_t::result_type;

    std::vector<hpx::future<result_t>> results;
    auto graph_partitions = partitions(G);
    results.reserve(graph_partitions.size());

    for (auto partition : graph_partitions) {
      results.push_back(dispatch_async(
        partition_locality(partition), algorithm_t{}, policy, std::true_type(),
        detail::remote_graph_ref(G),
        partition, ts...));
    }

    hpx::wait_all(results);
    return results;
  }

  /**
   * @brief Legacy HPX-segmented fallback for older distributed algorithms.
   *
   * This path relies on segmented iterator `.index()` semantics and is not part
   * of the generic partition-aware extension contract.
   */
  template <typename Algorithm, typename ExPolicy, detail::legacy_segmented_partitioned_graph Graph,
            typename... Ts>
    requires (!detail::is_partition_aware_algorithm<Algorithm>::value)
  [[gnu::noinline]] auto partitioned_segmented_algorithm(ExPolicy&& policy, Graph& G, Ts&&... ts) {
    using algorithm_t = Algorithm;
    using result_t = typename algorithm_t::result_type;

    auto first = G.begin();
    auto last = G.end();

    using traits = hpx::traits::segmented_iterator_traits<decltype(first.index())>;
    using segment_iterator = typename traits::segment_iterator;
    using local_iterator_type = typename traits::local_iterator;

    segment_iterator sit = traits::segment(first.index());
    segment_iterator send = traits::segment(last.index());

    std::vector<hpx::future<result_t>> results;
    results.reserve(send - sit + 1);

    if (sit == send) {
      // all elements are on the same partition
      local_iterator_type beg = traits::local(first.index());
      local_iterator_type end = traits::local(last.index());
      if (beg != end) {
        results.push_back(dispatch_async(
          traits::get_id(sit), algorithm_t{}, policy, std::true_type(),
          detail::remote_graph_ref(G),
          detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end), ts...));
      }
    }
    else {
      // handle all of or the remaining part of the first partition
      local_iterator_type beg = traits::local(first.index());
      local_iterator_type end = traits::end(sit);

      if (beg != end) {
        results.push_back(dispatch_async(
          traits::get_id(sit), algorithm_t{}, policy, std::true_type(),
          detail::remote_graph_ref(G),
          detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end), ts...));
      }

      // handle all full partitions except last
      for (++sit; sit != send; ++sit) {
        beg = traits::begin(sit);
        end = traits::end(sit);

        if (beg != end) {
          results.push_back(dispatch_async(
            traits::get_id(sit), algorithm_t{}, policy, std::true_type(),
            detail::remote_graph_ref(G),
            detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end), ts...));
        }
      }

      // handle the beginning (or all of) of the last partition
      beg = traits::begin(sit);
      end = traits::local(last.index());
      if (beg != end) {
        results.push_back(dispatch_async(
          traits::get_id(sit), algorithm_t{}, policy, std::true_type(),
          detail::remote_graph_ref(G),
          detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end), ts...));
      }
    }

    hpx::wait_all(results);
    return results;
  }
} // namespace nw::graph

#endif //  NW_GRAPH_PARTITIONED_ALGORITHM_HPP
