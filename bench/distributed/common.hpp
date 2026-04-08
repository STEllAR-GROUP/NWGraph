/**
 * @file common.hpp
 *
 * Distributed benchmark helpers.
 */

#ifndef NW_GRAPH_BENCH_DISTRIBUTED_COMMON_HPP
#define NW_GRAPH_BENCH_DISTRIBUTED_COMMON_HPP

#include "../common.hpp"

#include <utility>

#include <nwgraph/distributed/adjacency.hpp>
#include <nwgraph/distributed/serialize.hpp>

#include <hpx/include/partitioned_vector.hpp>

namespace nw::graph::bench {

inline auto load_partitioned_adjacency_graph(std::string mtx_file) {
  nw::util::life_timer _(__func__);
  std::string b_adj_file = partitioned_serialize_adj(mtx_file);
  return partitioned_deserialize_adj(b_adj_file);
}

using distributed_graph_type = decltype(load_partitioned_adjacency_graph(std::declval<std::string>()));
using distributed_degree_vector = hpx::partitioned_vector<typename distributed_graph_type::vertex_id_type>;
using distributed_rank_vector = hpx::partitioned_vector<float>;

} // namespace nw::graph::bench

#endif // NW_GRAPH_BENCH_DISTRIBUTED_COMMON_HPP