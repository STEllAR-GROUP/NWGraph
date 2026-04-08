#include <compare>
#include <cstdint>
#include <vector>

#include "../common/test_header.hpp"

#include <hpx/include/partitioned_vector.hpp>

#include "nwgraph/distributed/adjacency.hpp"

using registered_neighbors = std::vector<nw::graph::default_vertex_id_type>;
HPX_REGISTER_PARTITIONED_VECTOR(registered_neighbors)

using namespace nw::graph;

namespace usage_example {

using neighbor_list = std::vector<default_vertex_id_type>;

using graph_type = hpx::partitioned_vector<neighbor_list>;

} // namespace usage_example

namespace hpx {

hpx::id_type tag_invoke(nw::graph::vertex_partition_tag,
                        usage_example::graph_type const& G,
                        nw::graph::default_vertex_id_type v) {
  auto partition = G.get_segment_iterator(static_cast<std::size_t>(v));
  return hpx::naming::get_locality_from_id(partition->get_id());
}

nw::graph::default_vertex_id_type tag_invoke(nw::graph::target_tag,
                                             usage_example::graph_type const&,
                                             nw::graph::default_vertex_id_type v) {
  return v;
}

nw::graph::default_vertex_id_type tag_invoke(nw::graph::num_vertices_tag,
                                             usage_example::graph_type const& G) {
  return static_cast<nw::graph::default_vertex_id_type>(G.size());
}

} // namespace hpx

namespace nw::graph {

template <>
struct graph_traits<usage_example::graph_type> {
  using vertex_id_type = default_vertex_id_type;
};

} // namespace nw::graph

#include "nwgraph/distributed/algorithms/util.hpp"

namespace {

usage_example::graph_type make_partitioned_neighbor_graph() {
  std::vector<hpx::id_type> localities(2, hpx::find_here());
  usage_example::graph_type G(
    4u, hpx::explicit_container_layout(std::vector<std::size_t>{2u, 2u}, localities));

  G[0] = usage_example::neighbor_list{1u, 2u};
  G[1] = usage_example::neighbor_list{0u, 2u};
  G[2] = usage_example::neighbor_list{0u, 1u, 3u};
  G[3] = usage_example::neighbor_list{2u};

  return G;
}

template <typename Graph>
std::vector<vertex_id_t<std::remove_reference_t<Graph>>> targets_of(
  Graph& G, vertex_id_t<std::remove_reference_t<Graph>> u) {
  using graph_type = std::remove_reference_t<Graph>;
  using neighborhood_type = inner_range_t<graph_type>;
  std::vector<vertex_id_t<std::remove_reference_t<Graph>>> result;
  neighborhood_type neighborhood = G[u];
  for (auto&& edge : neighborhood) {
    result.push_back(target(G, edge));
  }
  return result;
}

} // namespace

TEST_CASE("partitioned algorithms accept a partitioned_vector range-of-ranges graph",
          "[distributed][partitioned][concepts][smoke]") {
  static_assert(graph<usage_example::graph_type>);
  static_assert(partitioned_graph<usage_example::graph_type>);
  static_assert(partitioned_algorithm_graph<usage_example::graph_type>);
  static_assert(std::same_as<partition_t<usage_example::graph_type>, hpx::id_type>);

  auto G = make_partitioned_neighbor_graph();

  // This example starts from hpx::partitioned_vector<std::vector<vertex>>.
  // The only distributed customization point it provides is vertex_partition(G, v),
  // which yields the HPX dispatch target directly.
  REQUIRE(vertex_partition(G, 0u) == hpx::find_here());
  REQUIRE(vertex_partition(G, 3u) == hpx::find_here());
  auto segments = partition_segments(G);
  REQUIRE(segments.size() == 2u);
  REQUIRE(segments[0].first == 0u);
  REQUIRE(segments[0].last == 2u);
  REQUIRE(segments[1].first == 2u);
  REQUIRE(segments[1].last == 4u);
  REQUIRE(is_local(G, 0u));
  REQUIRE(is_local(G, 2u));
  REQUIRE(targets_of(G, 0u) == (std::vector<default_vertex_id_type>{1u, 2u}));
  REQUIRE(targets_of(G, 2u) == (std::vector<default_vertex_id_type>{0u, 1u, 3u}));

  auto avg_degree = partitioned_avg_degree_per_partition(G);
  REQUIRE(avg_degree == (std::vector<double>{2.0, 2.0}));

  auto avg_remote_degree = partitioned_avg_remote_degree_per_partition(G);
  REQUIRE(avg_remote_degree == (std::vector<double>{1.0, 1.0}));
}
