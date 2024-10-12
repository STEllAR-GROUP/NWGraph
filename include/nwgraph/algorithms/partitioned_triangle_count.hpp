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

#ifndef NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_HPP
#define NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

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

    struct triangle_count : hpx::parallel::detail::algorithm<triangle_count, std::size_t> {

      constexpr triangle_count() noexcept
        : hpx::parallel::detail::algorithm<triangle_count, std::size_t>("triangle_count") {}

      template <typename ExPolicy, typename Graph>
      static std::size_t sequential(ExPolicy&& policy, Graph const& G, std::size_t first_index,
                                    std::size_t last_index) {
        std::size_t triangles = 0;
        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;
        for (auto u_neighbors = first; u_neighbors != last; ++u_neighbors) {
          for (auto elt = (*u_neighbors).begin(); elt != (*u_neighbors).end(); ++elt) {
            auto v = target(G, *elt);
            triangles += nw::graph::intersection_size(*u_neighbors, G[v]);
          }
        }
        return triangles;
      }

      template <typename ExPolicy, typename Graph, typename IterB, typename IterE>
      static std::size_t parallel(ExPolicy&& policy, Graph const& G, IterB first, IterE last) {
        return 0;
      }
    };
    /// \endcond

    template <typename Traits, typename SegIter, typename LocalIter>
    auto global_index(SegIter const& seg, LocalIter const& local) {

      auto composed = Traits::compose(seg, local);
      return composed.get_global_index();
    }
  } // namespace detail

  /**
   * @brief Two-dimensional triangle counting, parallel version.
   *
   * @tparam Graph adjacency_list_graph
   * @param G graph
   * @param threads number of threads
   * @return std::size_t number of triangles
   */
  template <typename ExPolicy, adjacency_list_graph Graph>
  [[gnu::noinline]] std::size_t partitioned_triangle_count(ExPolicy&& policy, Graph const& G) {
    auto first = G.begin();
    auto last = G.end();

    using traits = hpx::traits::segmented_iterator_traits<decltype(first.index())>;
    using segment_iterator = typename traits::segment_iterator;
    using local_iterator_type = typename traits::local_iterator;

    segment_iterator sit = traits::segment(first.index());
    segment_iterator send = traits::segment(last.index());

    using algorithm_t = detail::triangle_count;

    std::vector<hpx::future<std::size_t>> counts;
    counts.reserve(send - sit + 1);

    if (sit == send) {
      // all elements are on the same partition
      local_iterator_type beg = traits::local(first.index());
      local_iterator_type end = traits::local(last.index());
      if (beg != end) {
        counts.push_back(dispatch_async(
          traits::get_id(sit), algorithm_t{}, policy, std::true_type(), hpx::ref(G),
          detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end)));
      }
    }
    else {
      // handle all of or the remaining part of the first partition
      local_iterator_type beg = traits::local(first.index());
      local_iterator_type end = traits::end(sit);

      if (beg != end) {
        counts.push_back(dispatch_async(
          traits::get_id(sit), algorithm_t{}, policy, std::true_type(), hpx::ref(G),
          detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end)));
      }

      // handle all full partitions except last
      for (++sit; sit != send; ++sit) {
        beg = traits::begin(sit);
        end = traits::end(sit);

        if (beg != end) {
          counts.push_back(dispatch_async(
            traits::get_id(sit), algorithm_t{}, policy, std::true_type(), hpx::ref(G),
            detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end)));
        }
      }

      // handle the beginning (or all of) of the last partition
      beg = traits::begin(sit);
      end = traits::local(last.index());
      if (beg != end) {
        counts.push_back(dispatch_async(
          traits::get_id(sit), algorithm_t{}, policy, std::true_type(), hpx::ref(G),
          detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end)));
      }
    }

    hpx::wait_all(counts);
    return std::reduce(counts.begin(), counts.end(), std::size_t(0),
                       [](std::size_t count, auto&& f) { return count + f.get(); });
  }

  template <adjacency_list_graph Graph>
  std::size_t partitioned_triangle_count(Graph const& G) {
    return partitioned_triangle_count(hpx::execution::seq, G);
  }
} // namespace nw::graph

#endif //  NW_GRAPH_PARTITIONED_TRIANGLE_COUNT_HPP
