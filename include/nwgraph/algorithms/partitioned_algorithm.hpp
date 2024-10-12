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

#ifndef NW_GRAPH_PARTITIONED_ALGORITHM_HPP
#define NW_GRAPH_PARTITIONED_ALGORITHM_HPP

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

    template <typename Graph, typename Vertex>
    hpx::id_type vertex_locality(Graph const& G, Vertex v) {
      using traits = hpx::traits::segmented_iterator_traits<decltype(G.begin().index())>;
      return traits::get_id(traits::segment((G.begin() + v).index()));
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
  } // namespace detail

  /**
   * @brief Generic segemented algorithm implementation
   *
   * @tparam Graph adjacency_list_graph
   * @param G graph
   */
  template <typename Algorithm, typename ExPolicy, adjacency_list_graph Graph>
  [[gnu::noinline]] auto partitioned_algorithm(ExPolicy&& policy, Graph& G) {
    auto first = G.begin();
    auto last = G.end();

    using traits = hpx::traits::segmented_iterator_traits<decltype(first.index())>;
    using segment_iterator = typename traits::segment_iterator;
    using local_iterator_type = typename traits::local_iterator;

    segment_iterator sit = traits::segment(first.index());
    segment_iterator send = traits::segment(last.index());

    using algorithm_t = Algorithm;
    using result_t = typename algorithm_t::result_type;

    std::vector<hpx::future<result_t>> results;
    results.reserve(send - sit + 1);

    if (sit == send) {
      // all elements are on the same partition
      local_iterator_type beg = traits::local(first.index());
      local_iterator_type end = traits::local(last.index());
      if (beg != end) {
        results.push_back(dispatch_async(
          traits::get_id(sit), algorithm_t{}, policy, std::true_type(), hpx::ref(G),
          detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end)));
      }
    }
    else {
      // handle all of or the remaining part of the first partition
      local_iterator_type beg = traits::local(first.index());
      local_iterator_type end = traits::end(sit);

      if (beg != end) {
        results.push_back(dispatch_async(
          traits::get_id(sit), algorithm_t{}, policy, std::true_type(), hpx::ref(G),
          detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end)));
      }

      // handle all full partitions except last
      for (++sit; sit != send; ++sit) {
        beg = traits::begin(sit);
        end = traits::end(sit);

        if (beg != end) {
          results.push_back(dispatch_async(
            traits::get_id(sit), algorithm_t{}, policy, std::true_type(), hpx::ref(G),
            detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end)));
        }
      }

      // handle the beginning (or all of) of the last partition
      beg = traits::begin(sit);
      end = traits::local(last.index());
      if (beg != end) {
        results.push_back(dispatch_async(
          traits::get_id(sit), algorithm_t{}, policy, std::true_type(), hpx::ref(G),
          detail::global_index<traits>(sit, beg), detail::global_index<traits>(sit, end)));
      }
    }

    hpx::wait_all(results);
    return results;
  }
} // namespace nw::graph

#endif //  NW_GRAPH_PARTITIONED_ALGORITHM_HPP
