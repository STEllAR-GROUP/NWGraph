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

#include <nwgraph/containers/partitioned_compressed_local_view.hpp>
#include <nwgraph/partitioned_adjacency_local_view.hpp>
#include <nwgraph/util/partitioned_vector_local_partition_view.hpp>
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

    template <typename T>
    using local_pv_view = nw::graph::util::partitioned_vector_local_partition_view<T>;

    template <int idx, typename index_type, typename vertex_id, typename... Attributes>
    using local_adj_view =
      nw::graph::partitioned_index_adjacency_local_view<idx, index_type, vertex_id, Attributes...>;

    ////////////////////////////////////////////////////////////////////////////
    // handle counting of triangles on target locality
    template <typename Graph, typename vertex_id_t = typename Graph::vertex_id_type>
    static size_t triangle_counter_4(
      Graph G, size_t part_num,
      std::vector<std::tuple<std::vector<vertex_id_t>, std::vector<vertex_id_t>>> const&
        packets) {

      hpx::scoped_annotation ann_tc_handle_packet("TC_handle_packet");

      local_adj_view G_loc(G, part_num);

      auto get_tgt = [&G_loc](auto&& e) { return target(G_loc, e); };

      size_t triangles = 0;
      for (auto& packet : packets) {
        for (auto& v : std::get<0>(packet)) {
          auto& set1 = std::get<1>(packet);
          auto set2 = G_loc[v] | std::ranges::views ::transform(get_tgt);
          triangles +=
            nw::graph::intersection_size(set1, set2, std::less<vertex_id_t>{});
        }
      }
      return triangles;
    }

    template <typename Graph>
    struct triangle_count_action_4
      : hpx::actions::action<decltype(&triangle_counter_4<Graph>), &triangle_counter_4<Graph>,
                             triangle_count_action_4<Graph>> {};

    ////////////////////////////////////////////////////////////////////////////
    struct triangle_count_4 : hpx::parallel::detail::algorithm<triangle_count_4, size_t> {

      // triangle counting driver for one of the partitions
      constexpr triangle_count_4() noexcept
        : hpx::parallel::detail::algorithm<triangle_count_4, size_t>("triangle_count_4") {}


      template <typename ExPolicy, typename Graph>
      static size_t sequential(ExPolicy&&, Graph G, size_t first_index, size_t last_index,
                             size_t batchsize) {
        // Convert to local views and call the main implementation
        std::size_t partnum = G.get_indices().get_partition(first_index);
        local_adj_view G_loc(G, partnum);
        return seq_impl(G_loc, batchsize);
      }

      template <typename LocGraph>
      static size_t seq_impl(LocGraph G_loc, size_t batchsize) {

        hpx::scoped_annotation ann_tc_impl("TC_impl");

        using Graph = LocGraph::graph_type;
        using vertex_id_type = typename Graph::vertex_id_type;
        using target_list_t = std::vector<
          std::tuple<std::vector<vertex_id_type>, std::vector<vertex_id_type>>>;
        using remote_targets_t = std::map<std::size_t, target_list_t>;

        auto cmp = [&G_loc](auto&& a, auto&& b) { return target(G_loc, a) < target(G_loc, b); };

        auto get_part_locality = [](auto& graph_loc, std::size_t part_num) -> hpx::id_type
        {
          // TODO: Too intrusive, fix
          hpx::id_type part_id = graph_loc.parent().get_indices().partitions()[part_num].get_id();
          return hpx::naming::get_locality_from_id(part_id);
        };

        auto send_remote_action = [&](std::size_t part_num, auto&& targets) -> hpx::future<size_t>
        {
          hpx::scoped_annotation ann_tc_send_remote("TC_send_remote");
          triangle_count_action_4<Graph> act;
          auto id = get_part_locality(G_loc, part_num);
          return hpx::async(act, id, hpx::ref(G_loc.parent()), part_num, std::move(targets));
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
          std::map<std::size_t, std::vector<vertex_id_type>> target_vertices;

          for (auto const& edge : neighbor_range) {

            vertex_id_type v = target(G_loc, edge);

            if (G_loc.is_local_index(v)) {
              // handle things locally
              triangles += nw::graph::intersection_size(neighbor_range, G_loc[v], cmp);
            }
            else {
              // send our neighbors to elt's locality
              auto part_num = vertex_partition_num(G_loc.parent(), v);
              target_vertices[part_num].push_back(v);
            }

            // collect all neighbor vertex ids for current vertex
            neighbors.push_back(v);
          }

          counts.get() += triangles;

          // now store this vertex list in collection of messages to send
          for (auto&& [part_num, vertices] : target_vertices) {
            auto& targets = remote_targets.get()[part_num];
            // store list of vertices with their neighbors in any case
            targets.push_back(std::make_tuple(std::move(vertices), neighbors));

            // if batch size has been reached, trigger async operation
            if (targets.size() >= batchsize) {
              auto fut = send_remote_action(part_num, std::move(targets));
              targets = target_list_t{};
              remote_results.get().push_back(std::move(fut));
            }
          }
          
        };

        {
          hpx::scoped_annotation ann_tc_main_loop("TC_main_loop");
          hpx::for_each(hpx::execution::par, G_loc.begin(), G_loc.end(), tc);
        }

        // Send any remaining messages
        remote_targets.reduce(
          [&](auto&& thd_remote_targets)
          {
            for (auto& [part_num, targets] : thd_remote_targets) {
              if (!targets.empty()) {
                auto fut = send_remote_action(part_num, std::move(targets));
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

  template <adjacency_list_graph Graph>
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
