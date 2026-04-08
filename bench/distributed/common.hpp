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
#include <nwgraph/distributed/algorithms/util.hpp>
#include <nwgraph/distributed/serialize.hpp>

#include <hpx/include/partitioned_vector.hpp>
#include <hpx/include/partitioned_vector_predef.hpp>

namespace nw::graph::bench {

template <int idx, directedness dir>
inline auto load_partitioned_adjacency_graph(std::string mtx_file) {
  nw::util::life_timer _(__func__);
  std::string b_adj_file = partitioned_serialize_adj<idx, dir>(mtx_file);
  return partitioned_deserialize_adj<idx, dir>(b_adj_file);
}

template <int idx, directedness dir>
using distributed_graph_type = decltype(load_partitioned_adjacency_graph<idx, dir>(std::declval<std::string>()));

using distributed_degree_vector = hpx::partitioned_vector<typename distributed_graph_type<0, directedness::directed>::vertex_id_type>;
using distributed_rank_vector = hpx::partitioned_vector<float>;

template <partitioned_algorithm_graph Graph>
auto build_random_sources(Graph& graph, size_t n, long seed) {
  using Id = typename nw::graph::vertex_id_t<std::remove_reference_t<Graph>>;
  using traits = hpx::traits::segmented_iterator_traits<decltype(std::declval<hpx::partitioned_vector<Id>&>().begin())>;

  auto sources = std::vector<Id>(n);
  auto degrees = partitioned_row_degrees(graph);

  std::vector<Id> local_degrees;
  local_degrees.reserve(degrees.size());

  std::size_t num_partitions = traits::segment(degrees.end()) - traits::segment(degrees.begin());
  for (std::size_t part = 0; part != num_partitions; ++part) {
    auto values = degrees.get_values(hpx::launch::sync, part);
    std::move(values.begin(), values.end(), std::back_inserter(local_degrees));
  }

  auto gen = std::mt19937(seed);
  auto dis = std::uniform_int_distribution<Id>(0, num_vertices(graph) - 1);

  for (auto& id : sources) {
    for (id = dis(gen); local_degrees[id] == 0; id = dis(gen)) {
    }
  }
  return sources;
}

} // namespace nw::graph::bench

#endif // NW_GRAPH_BENCH_DISTRIBUTED_COMMON_HPP