/**
 * @file partitioned_bfs_1.hpp
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

#ifndef NW_GRAPH_PARTITIONED_BFS_1_HPP
#define NW_GRAPH_PARTITIONED_BFS_1_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/algorithms/bfs.hpp"
#include "nwgraph/distributed/algorithms/algorithm.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <vector>

#include <hpx/executors/execution_policy.hpp>
#include <hpx/include/partitioned_vector_predef.hpp>
#include <hpx/parallel/util/detail/algorithm_result.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace nw::graph {

  namespace detail {

    ////////////////////////////////////////////////////////////////////////////
    template <typename T>
    auto copy_to_local(hpx::partitioned_vector<T> const& d) {

      using traits = hpx::traits::segmented_iterator_traits<decltype(d.begin())>;

      std::vector<T> result;
      result.reserve(d.size());

      // handle partitions separately
      std::size_t num_partitions = traits::segment(d.end()) - traits::segment(d.begin());
      for (std::size_t part = 0; part != num_partitions; ++part) {
        auto data = d.get_values(hpx::launch::sync, part);
        std::move(data.begin(), data.end(), std::back_inserter(result));
      }

      return result;
    }

    ////////////////////////////////////////////////////////////////////////////
    template <typename Vertex>
    bool set_parent(hpx::partitioned_vector<Vertex>& parents, Vertex v, Vertex u, size_t) {
      using traits = hpx::traits::segmented_iterator_traits<decltype(parents.begin())>;

      // *traits::local() returns a proxy object, local_get() returns the
      // reference to the element (assuming it is local)
      std::atomic_ref<Vertex> parent((*traits::local(parents.begin() + v)).local_get());
      Vertex curr_parent = parent.load(std::memory_order_acquire);

      if (curr_parent == std::numeric_limits<Vertex>::max()) {
        while (!parent.compare_exchange_strong(curr_parent, u,
                                               std::memory_order::acq_rel)) {
          if (u >= curr_parent) {
            // since we last checked, some other thread has set a better parent
            return false;
          }
        }
        return true; // u is now the new parent of v, keep traversing v
      }

      // v already has a parent
      return false;
    }

    ////////////////////////////////////////////////////////////////////////////
    // handle counting of triangles on target locality
    template <typename Vertex>
    struct bfs_action_1;

    template <typename Graph>
    static void
    bfs_1(Graph G, hpx::partitioned_vector<vertex_id_t<Graph>> parents, size_t batchsize,
          std::vector<std::tuple<vertex_id_t<Graph>, vertex_id_t<Graph>, size_t>> const& sources) {

      std::uint32_t this_locality_id = hpx::get_locality_id();

      using vertex_id_type = vertex_id_t<Graph>;
      using target_list_t = std::vector<std::tuple<vertex_id_type, vertex_id_type, size_t>>;
      using remote_traversal_t = std::map<std::uint32_t, target_list_t>;

      // handle source vertices separately
      std::deque<std::tuple<vertex_id_type, size_t>> q1, q2;
      for (auto const& source : sources) {
        assert(is_same_locality(this_locality_id, G, std::get<0>(source)));
        if (set_parent(parents, std::get<0>(source), std::get<1>(source), std::get<2>(source))) {
          // traverse into neighbors only if source was not traversed before
          q1.push_back(std::make_tuple(std::get<0>(source), std::get<2>(source)));
        }
      }

      // for each v in G do
      remote_traversal_t remote_traversals;
      std::vector<hpx::future<void>> remote_ops;

      auto tc = [&](std::tuple<vertex_id_type, size_t> const& next)
      {
        auto [u, lvl] = next;
        for (auto const& edge : G[u]) {
          vertex_id_type v = target(G, edge);
          auto id = vertex_locality_id(G, v);

          if (id == this_locality_id) {
            // handle things locally
            if (set_parent(parents, v, u, lvl)) {
              // traverse into neighbors only if v was not traversed before
              q2.push_back(std::make_tuple(v, lvl + 1));
            }
          }
          else {
            // send traversal request to correct locality
            auto& targets = remote_traversals[id];

            // if batch size has been reached, trigger async operation
            if (targets.size() >= batchsize) {
              bfs_action_1<Graph> act;
              remote_ops.push_back(hpx::async(act, hpx::naming::get_id_from_locality_id(id),
                                              hpx::ref(G), hpx::ref(parents), batchsize,
                                              std::move(targets)));
              targets = target_list_t{};
            }

            // store list of vertices to traverse remotely in any case
            targets.push_back(std::make_tuple(v, u, lvl + 1));
          }
        }
      };

      // do actual traversal
      while (!q1.empty()) {
        std::for_each(q1.begin(), q1.end(), tc);
        std::swap(q1, q2);
        q2.clear();
      }

      // trigger remaining asynchronous operations
      for (auto&& [id, targets] : remote_traversals) {
        if (!targets.empty()) {
          bfs_action_1<Graph> act;
          remote_ops.push_back(hpx::async(act, hpx::naming::get_id_from_locality_id(id),
                                          hpx::ref(G), hpx::ref(parents), batchsize,
                                          std::move(targets)));
        }
      }

      // wait for all remote operations to finish
      if (!remote_ops.empty()) {
        hpx::wait_all(remote_ops);
      }
    }

    template <typename Graph>
    struct bfs_action_1
      : hpx::actions::action<decltype(&bfs_1<Graph>), &bfs_1<Graph>, bfs_action_1<Graph>> {};
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
  auto partitioned_bfs_1(Graph& G, vertex_id_t<Graph> source, size_t batchsize) {

    auto sizes = G.indices_.get_partition_sizes();
    size_t size = std::reduce(sizes.begin(), sizes.end(), size_t(0)) - 1;

    using vertex_id_type = vertex_id_t<Graph>;

    hpx::partitioned_vector<vertex_id_type> parents(
      size, std::numeric_limits<vertex_id_type>::max(),
      hpx::explicit_container_layout(std::move(sizes), G.indices_.get_partition_localities()));
    parents.register_as("parents");

    std::vector<std::tuple<vertex_id_type, vertex_id_type, size_t>> targets;
    targets.push_back(std::make_tuple(source, source, 1));

    detail::bfs_action_1<Graph> act;
    auto f = hpx::async(act, detail::vertex_locality(G, source), hpx::ref(G), hpx::ref(parents),
                        batchsize, std::move(targets));

    hpx::wait_all(f);

    return detail::copy_to_local(parents);
  }
} // namespace nw::graph

#endif //  NW_GRAPH_PARTITIONED_BFS_1_HPP
