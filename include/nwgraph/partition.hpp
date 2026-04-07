/**
 * @file partition.hpp
 *
 * Partition-centric helpers for distributed NWGraph containers.
 */

#ifndef NW_GRAPH_PARTITION_HPP
#define NW_GRAPH_PARTITION_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/partitioned_adjacency.hpp"
#include "nwgraph/partitioned_adjacency_local_view.hpp"
#include "nwgraph/util/partitioned_vector_local_partition_view.hpp"

#include <cstddef>
#include <type_traits>

#include <hpx/assert.hpp>
#include <hpx/include/partitioned_vector.hpp>

namespace nw::graph {

  template <typename T>
  auto tag_invoke(local_view_tag, hpx::partitioned_vector<T>& pv, partition_descriptor partition) {
    return util::partitioned_vector_local_partition_view<T>(pv, partition);
  }

  template <typename T>
  auto tag_invoke(remote_ref_tag, hpx::partitioned_vector<T>& pv) {
    using partitions_type = std::decay_t<decltype(pv.partitions())>;
    return hpx::partitioned_vector<T>::create_from(
      hpx::id_type{}, pv.size(), partitions_type(pv.partitions()));
  }

  template <typename T, std::unsigned_integral Index>
  bool tag_invoke(is_local_index_tag, util::partitioned_vector_local_partition_view<T> const& pv,
                  Index i) {
    return pv.is_local_index(static_cast<std::size_t>(i));
  }

  template <typename T>
  partition_descriptor tag_invoke(partition_tag,
                                  util::partitioned_vector_local_partition_view<T> const& pv) {
    return pv.partition();
  }

  template <typename T, copartitioned_graph Graph>
  hpx::partitioned_vector<T> make_copartitioned_vector(Graph const& G, T init_value = T{}) {
    auto const& indices = G.get_indices();
    return hpx::partitioned_vector<T>(
      indices.size(), init_value,
      hpx::explicit_container_layout(indices.get_partition_sizes(), indices.get_partition_localities()));
  }

} // namespace nw::graph

#endif // NW_GRAPH_PARTITION_HPP
