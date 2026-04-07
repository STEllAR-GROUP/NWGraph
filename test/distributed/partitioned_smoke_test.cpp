#include <atomic>
#include <compare>
#include <cstdint>
#include <string>
#include <vector>

#include "../common/test_header.hpp"

#include <hpx/include/partitioned_vector.hpp>

#include "nwgraph/partition.hpp"

using registered_uint = std::uint32_t;
HPX_REGISTER_PARTITIONED_VECTOR(registered_uint)

using namespace nw::graph;

namespace {

std::string unique_graph_name(char const* prefix) {
  static std::atomic<std::uint32_t> counter{0};
  return std::string(prefix) + "_" + std::to_string(counter++);
}

partitioned_adjacency<0> make_partitioned_smoke_graph() {
  // Build a tiny local graph first, then distribute it with an explicit
  // vertex partitioning and edge partitioning. This is the simplest
  // end-to-end construction path for the partitioned graph type.
  edge_list<directedness::directed> edges{
    {0u, 1u},
    {0u, 2u},
    {1u, 0u},
    {1u, 2u},
    {2u, 0u},
    {2u, 1u},
    {2u, 3u},
    {3u, 2u},
  };

  adjacency<0> local_graph(edges);

  std::vector<std::size_t> vertex_partition_sizes{2, 2};
  std::vector<std::size_t> edge_partition_sizes{4, 4};
  std::vector<hpx::id_type> localities(2, hpx::find_here());
  auto name = unique_graph_name("partitioned_smoke");

  return partitioned_adjacency<0>(
    4, 8, vertex_partition_sizes, edge_partition_sizes, local_graph, name.c_str(), localities);
}

template <typename Graph>
std::vector<typename std::remove_reference_t<Graph>::vertex_id_type> targets_of(
  Graph& G, typename std::remove_reference_t<Graph>::vertex_id_type u) {
  std::vector<typename std::remove_reference_t<Graph>::vertex_id_type> result;
  for (auto&& edge : G[u]) {
    result.push_back(target(G, edge));
  }
  return result;
}

template <typename Index>
std::vector<Index> neighbors_of(
  hpx::partitioned_vector<Index>& offsets, hpx::partitioned_vector<Index>& neighbors,
  std::size_t u) {
  auto begin = static_cast<Index>(offsets[u]);
  auto end = static_cast<Index>(offsets[u + 1]);

  std::vector<Index> result;
  for (auto i = begin; i < end; ++i) {
    result.push_back(static_cast<Index>(neighbors[i]));
  }
  return result;
}

} // namespace

TEST_CASE("partitioned adjacency smoke test", "[distributed][partitioned][smoke]") {
  auto G = make_partitioned_smoke_graph();

  // Discover the logical graph partitions and ask which vertices they own.
  auto graph_partitions = partitions(G);
  REQUIRE(graph_partitions.size() == 2);

  REQUIRE(graph_partitions[0].first_index() == 0u);
  REQUIRE(graph_partitions[0].last_index() == 2u);
  REQUIRE(graph_partitions[1].first_index() == 2u);
  REQUIRE(graph_partitions[1].last_index() == 4u);

  REQUIRE(partition_locality(graph_partitions[0]) == hpx::find_here());
  REQUIRE(partition_locality(graph_partitions[1]) == hpx::find_here());

  REQUIRE(vertex_partition(G, 0u) == graph_partitions[0]);
  REQUIRE(vertex_partition(G, 1u) == graph_partitions[0]);
  REQUIRE(vertex_partition(G, 2u) == graph_partitions[1]);
  REQUIRE(vertex_partition(G, 3u) == graph_partitions[1]);

  REQUIRE(is_local(graph_partitions[0], 0u));
  REQUIRE(is_local(graph_partitions[0], 1u));
  REQUIRE_FALSE(is_local(graph_partitions[0], 2u));
  REQUIRE(is_local(graph_partitions[1], 2u));
  REQUIRE(is_local(graph_partitions[1], 3u));

  // The distributed graph itself still behaves like an adjacency-list graph.
  REQUIRE(targets_of(G, 0u) == (std::vector<default_vertex_id_type>{1u, 2u}));
  REQUIRE(targets_of(G, 1u) == (std::vector<default_vertex_id_type>{0u, 2u}));
  REQUIRE(targets_of(G, 2u) == (std::vector<default_vertex_id_type>{0u, 1u, 3u}));
  REQUIRE(targets_of(G, 3u) == (std::vector<default_vertex_id_type>{2u}));

  // Open local graph views using the semantic partition descriptors.
  auto G0 = local_view(G, graph_partitions[0]);
  auto G1 = local_view(G, graph_partitions[1]);

  REQUIRE(G0.size() == 2u);
  REQUIRE(G1.size() == 2u);
  REQUIRE(is_local_index(G0, 0u));
  REQUIRE(is_local_index(G0, 1u));
  REQUIRE_FALSE(is_local_index(G0, 2u));
  REQUIRE(is_local_index(G1, 2u));
  REQUIRE(is_local_index(G1, 3u));
  REQUIRE(targets_of(G0, 0u) == (std::vector<default_vertex_id_type>{1u, 2u}));
  REQUIRE(targets_of(G0, 1u) == (std::vector<default_vertex_id_type>{0u, 2u}));
  REQUIRE(targets_of(G1, 2u) == (std::vector<default_vertex_id_type>{0u, 1u, 3u}));
  REQUIRE(targets_of(G1, 3u) == (std::vector<default_vertex_id_type>{2u}));

  // A plain partitioned_vector can be made to follow the same partitioning,
  // which is how algorithms keep auxiliary data aligned with the graph.
  std::vector<hpx::id_type> localities(2, hpx::find_here());
  hpx::partitioned_vector<default_index_t> values(
    4u, 7u, hpx::explicit_container_layout(std::vector<std::size_t>{2u, 2u}, localities));
  auto values_0 = local_view(values, graph_partitions[0]);
  auto values_1 = local_view(values, graph_partitions[1]);

  REQUIRE(is_local_index(values_0, 0u));
  REQUIRE(is_local_index(values_0, 1u));
  REQUIRE_FALSE(is_local_index(values_0, 2u));
  REQUIRE(is_local_index(values_1, 2u));
  REQUIRE(is_local_index(values_1, 3u));

  values_0[0u] = 9u;
  values_0[1u] = 11u;
  values_1[2u] = 13u;
  values_1[3u] = 15u;

  REQUIRE(static_cast<default_index_t>(values[0]) == 9u);
  REQUIRE(static_cast<default_index_t>(values[1]) == 11u);
  REQUIRE(static_cast<default_index_t>(values[2]) == 13u);
  REQUIRE(static_cast<default_index_t>(values[3]) == 15u);
}

TEST_CASE("partitioned_vector range-of-ranges smoke test", "[distributed][partitioned_vector][smoke]") {
  // This is the low-level "range of ranges" view of a graph:
  // offsets[u]..offsets[u+1] names the neighbor range for vertex u.
  std::vector<hpx::id_type> localities(2, hpx::find_here());

  hpx::partitioned_vector<default_index_t> offsets(
    5u, hpx::explicit_container_layout(std::vector<std::size_t>{3u, 2u}, localities));
  hpx::partitioned_vector<default_index_t> neighbors(
    8u, hpx::explicit_container_layout(std::vector<std::size_t>{4u, 4u}, localities));

  std::vector<default_index_t> offsets_data{0u, 2u, 4u, 7u, 8u};
  std::vector<default_index_t> neighbors_data{1u, 2u, 0u, 2u, 0u, 1u, 3u, 2u};

  for (std::size_t i = 0; i < offsets_data.size(); ++i) {
    offsets[i] = offsets_data[i];
  }
  for (std::size_t i = 0; i < neighbors_data.size(); ++i) {
    neighbors[i] = neighbors_data[i];
  }

  // HPX exposes the outer partitioning directly through segment iterators.
  auto first_segment = offsets.get_segment_iterator(0u);
  auto second_segment = offsets.get_segment_iterator(3u);

  REQUIRE(static_cast<std::size_t>(first_segment->first_) == 0u);
  REQUIRE(static_cast<std::size_t>(first_segment->size_) == 3u);
  REQUIRE(static_cast<std::size_t>(second_segment->first_) == 3u);
  REQUIRE(static_cast<std::size_t>(second_segment->size_) == 2u);

  // Interpreting the two vectors together gives the same graph as above.
  REQUIRE(neighbors_of(offsets, neighbors, 0u) == (std::vector<default_index_t>{1u, 2u}));
  REQUIRE(neighbors_of(offsets, neighbors, 1u) == (std::vector<default_index_t>{0u, 2u}));
  REQUIRE(neighbors_of(offsets, neighbors, 2u) == (std::vector<default_index_t>{0u, 1u, 3u}));
  REQUIRE(neighbors_of(offsets, neighbors, 3u) == (std::vector<default_index_t>{2u}));
}
