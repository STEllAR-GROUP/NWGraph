/**
 * @file copartitioned_vector.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef NW_GRAPH_UTIL_COPARTITIONED_VECTOR_HPP
#define NW_GRAPH_UTIL_COPARTITIONED_VECTOR_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/distributed/graph.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include <hpx/include/partitioned_vector.hpp>

namespace nw::graph {

  template <typename T, partitioned_graph Graph>
  hpx::partitioned_vector<T> make_copartitioned_vector(Graph const& G, T init_value = T{}) {
    auto segments = partition_segments(G);
    std::vector<std::size_t> sizes;
    std::vector<hpx::id_type> localities;
    sizes.reserve(segments.size());
    localities.reserve(segments.size());

    for (auto const& segment : segments) {
      sizes.push_back(segment.last - segment.first);
      localities.push_back(segment.locality);
    }

    return hpx::partitioned_vector<T>(
      static_cast<std::size_t>(G.size()), init_value,
      hpx::explicit_container_layout(std::move(sizes), std::move(localities)));
  }

} // namespace nw::graph

#endif // NW_GRAPH_UTIL_COPARTITIONED_VECTOR_HPP