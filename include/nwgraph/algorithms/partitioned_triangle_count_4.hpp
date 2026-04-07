/**
 * @file partitioned_triangle_count_4.hpp
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

#ifndef NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_4_HPP
#define NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_4_HPP

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
    // handle counting of triangles on target locality
    template <typename Graph, typename Partition>
    static size_t triangle_counter_4(
      Graph G, Partition partition,
      std::vector<std::tuple<std::vector<vertex_id_t<std::remove_reference_t<Graph>>>,
                             std::vector<vertex_id_t<std::remove_reference_t<Graph>>>>> const&
        packets) {

      hpx::scoped_annotation ann_tc_handle_packet("TC_handle_packet");

      decltype(auto) G_data = detail::partition_data(G, partition);

      auto get_tgt = [&G_data](auto&& e) { return target(G_data, e); };

      size_t triangles = 0;
      for (auto& packet : packets) {
        for (auto& v : std::get<0>(packet)) {
          auto& set1 = std::get<1>(packet);
          auto set2 = G_data[v] | std::ranges::views ::transform(get_tgt);
          triangles +=
            nw::graph::intersection_size(
              set1, set2, std::less<vertex_id_t<std::remove_reference_t<Graph>>>{});
        }
      }
      return triangles;
    }

    template <typename Graph, typename Partition>
    struct triangle_count_action_4
      : hpx::actions::action<decltype(&triangle_counter_4<Graph, Partition>),
                             &triangle_counter_4<Graph, Partition>,
                             triangle_count_action_4<Graph, Partition>> {};

    ////////////////////////////////////////////////////////////////////////////
    struct triangle_count_4 : hpx::parallel::detail::algorithm<triangle_count_4, size_t> {
      static constexpr bool partition_aware = true;

      // triangle counting driver for one of the partitions
      constexpr triangle_count_4() noexcept
        : hpx::parallel::detail::algorithm<triangle_count_4, size_t>("triangle_count_4") {}


      template <typename ExPolicy, typename Graph>
      static size_t sequential(ExPolicy&&, Graph G, partition_t<std::remove_reference_t<Graph>> partition,
                             size_t batchsize) {
        decltype(auto) G_data = detail::partition_data(G, partition);

        hpx::scoped_annotation ann_tc_impl("TC_impl");

        using graph_type = std::remove_reference_t<Graph>;
        using vertex_id_type = vertex_id_t<graph_type>;
        using partition_type = partition_t<graph_type>;
        using target_list_t = std::vector<
          std::tuple<std::vector<vertex_id_type>, std::vector<vertex_id_type>>>;
        using remote_targets_t = std::map<partition_type, target_list_t>;

        auto cmp = [&G_data](auto&& a, auto&& b) { return target(G_data, a) < target(G_data, b); };

        auto send_remote_action = [&](partition_type const& target_partition,
                                      auto&& targets) -> hpx::future<size_t>
        {
          hpx::scoped_annotation ann_tc_send_remote("TC_send_remote");
          triangle_count_action_4<graph_type, partition_type> act;
          return hpx::async(
            act, partition_locality(target_partition), detail::remote_graph_ref(G),
            target_partition, std::move(targets));
        };

        // for each v in G do
        safe_object<size_t> counts;
        safe_object<remote_targets_t> remote_targets;
        safe_object<std::vector<hpx::future<size_t>>> remote_results;

        auto tc = [&](auto&& neighbor_range)
        {
          hpx::scoped_annotation ann_tc_per_vertex("TC_per_vertex");
          size_t triangles = 0;
          std::vector<vertex_id_type> neighbors;
          std::map<partition_type, std::vector<vertex_id_type>> target_vertices;

          for (auto const& edge : neighbor_range) {

            vertex_id_type v = detail::edge_target(G_data, edge);

            if (is_local(partition, v)) {
              // handle things locally
              auto target_neighbors = detail::iterable_row(G_data[v]);
              triangles += nw::graph::intersection_size(neighbor_range, target_neighbors, cmp);
            }
            else {
              // send our neighbors to elt's locality
              auto target_partition = vertex_partition(G, v);
              target_vertices[target_partition].push_back(v);
            }

            // collect all neighbor vertex ids for current vertex
            neighbors.push_back(v);
          }

          counts.get() += triangles;

          // now store this vertex list in collection of messages to send
          for (auto&& [target_partition, vertices] : target_vertices) {
            auto& targets = remote_targets.get()[target_partition];
            // store list of vertices with their neighbors in any case
            targets.push_back(std::make_tuple(std::move(vertices), neighbors));

            // if batch size has been reached, trigger async operation
            if (targets.size() >= batchsize) {
              auto fut = send_remote_action(target_partition, std::move(targets));
              targets = target_list_t{};
              remote_results.get().push_back(std::move(fut));
            }
          }
          
        };

        {
          hpx::scoped_annotation ann_tc_main_loop("TC_main_loop");
          hpx::experimental::for_loop(
            hpx::execution::par,
            static_cast<vertex_id_type>(partition_first_index(partition)),
            static_cast<vertex_id_type>(partition_last_index(partition)),
            [&](vertex_id_type u) {
              auto neighbor_range = detail::iterable_row(G_data[u]);
              tc(neighbor_range);
            });
        }

        // Send any remaining messages
        remote_targets.reduce(
          [&](auto&& thd_remote_targets)
          {
            for (auto& [target_partition, targets] : thd_remote_targets) {
              if (!targets.empty()) {
                auto fut = send_remote_action(target_partition, std::move(targets));
                remote_results.get().push_back(std::move(fut));
              }
            }
          });


        size_t triangles = 0;

        // Accumulate local results
        counts.reduce([&](auto&& thd_counts) { triangles += thd_counts; });


        // Now wait/accumulate all remote results
        {
          hpx::scoped_annotation ann_tc_wait_for_results("TC_wait_for_results");
          remote_results.reduce(
            [&](auto&& thd_remote_results)
            {
              for (auto&& f : thd_remote_results) {
                triangles += f.get();
              }
            });
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

  template <partitioned_algorithm_graph Graph>
  size_t partitioned_triangle_count_4(Graph& G, size_t batchsize) {
    hpx::scoped_annotation ann_tc("TC");
    auto counts =
      partitioned_algorithm<detail::triangle_count_4>(hpx::execution::seq, G, batchsize);
    return std::transform_reduce(
      counts.begin(), counts.end(), size_t(0),
      [](size_t count, size_t curr) { return count + curr; }, [](auto&& f) { return f.get(); });
  }
} // namespace nw::graph

#endif //  NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_4_HPP
