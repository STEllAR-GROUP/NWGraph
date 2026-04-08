/**
 * @file partitioned_graph.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef NW_GRAPH_PARTITIONED_GRAPH_HPP
#define NW_GRAPH_PARTITIONED_GRAPH_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/graph_concepts.hpp"

#include <concepts>
#include <cstddef>
#include <ranges>
#include <type_traits>

#include <hpx/runtime_distributed/find_here.hpp>

namespace nw::graph {

  struct partition_segment {
    hpx::id_type locality;
    std::size_t first;
    std::size_t last;
  };

  DECL_TAG_INVOKE(vertex_partition);
  DECL_TAG_INVOKE(partition_segments);

  namespace detail {
    template <typename Segments>
    concept partition_segment_range =
      std::ranges::input_range<Segments> &&
      requires(std::ranges::range_reference_t<Segments> segment) {
        { segment.locality } -> std::convertible_to<hpx::id_type>;
        { segment.first } -> std::convertible_to<std::size_t>;
        { segment.last } -> std::convertible_to<std::size_t>;
      };
  } // namespace detail

  template <typename G>
  using partition_t = hpx::id_type;

  template <typename Graph, std::unsigned_integral VertexId>
  bool is_local(Graph const& G, VertexId v) {
    return vertex_partition(G, v) == hpx::find_here();
  }

  template <typename G>
  concept partitioned_graph =
    requires(std::remove_reference_t<G> const& g,
             vertex_id_t<std::remove_reference_t<G>> v) {
      typename vertex_id_t<std::remove_reference_t<G>>;
      { g.size() } -> std::convertible_to<std::size_t>;
      { vertex_partition(g, v) } -> std::same_as<hpx::id_type>;
      requires detail::partition_segment_range<decltype(partition_segments(g))>;
    };

} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_GRAPH_HPP