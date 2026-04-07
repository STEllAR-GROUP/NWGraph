#include <compare>
#include <cstdint>
#include <vector>

#include "../common/test_header.hpp"

#include <hpx/include/partitioned_vector.hpp>

#include "nwgraph/partition.hpp"

using registered_neighbors = std::vector<nw::graph::default_vertex_id_type>;
HPX_REGISTER_PARTITIONED_VECTOR(registered_neighbors)

using namespace nw::graph;

namespace usage_example {

using neighbor_list = std::vector<default_vertex_id_type>;

using graph_type = hpx::partitioned_vector<neighbor_list>;

struct simple_partition {
  hpx::id_type locality;
  std::size_t first = 0;
  std::size_t last = 0;

  auto operator<=>(simple_partition const&) const = default;

  template <typename Archive>
  void serialize(Archive& ar, unsigned) {
    ar& locality& first& last;
  }
};

inline hpx::id_type tag_invoke(nw::graph::partition_locality_tag,
                               simple_partition const& partition) {
  return partition.locality;
}

inline std::size_t tag_invoke(nw::graph::partition_first_index_tag,
                              simple_partition const& partition) {
  return partition.first;
}

inline std::size_t tag_invoke(nw::graph::partition_last_index_tag,
                              simple_partition const& partition) {
  return partition.last;
}

} // namespace usage_example

namespace hpx {

auto tag_invoke(nw::graph::partitions_tag, usage_example::graph_type const& G) {
  std::vector<usage_example::simple_partition> result;
  auto partition = G.segment_begin();
  auto partition_end = G.segment_end();
  for (; partition != partition_end; ++partition) {
    auto locality = hpx::naming::get_locality_from_id(partition->get_id());
    auto first_index = static_cast<std::size_t>(partition->first_);
    auto last_index = static_cast<std::size_t>(partition->first_ + partition->size_);
    last_index = std::min(last_index, static_cast<std::size_t>(G.size()));
    result.push_back({HPX_MOVE(locality), first_index, last_index});
  }
  return result;
}

usage_example::simple_partition tag_invoke(nw::graph::vertex_partition_tag,
                                           usage_example::graph_type const& G,
                                           nw::graph::default_vertex_id_type v) {
  auto partition = G.get_segment_iterator(static_cast<std::size_t>(v));
  auto locality = hpx::naming::get_locality_from_id(partition->get_id());
  auto first_index = static_cast<std::size_t>(partition->first_);
  auto last_index = static_cast<std::size_t>(partition->first_ + partition->size_);
  last_index = std::min(last_index, static_cast<std::size_t>(G.size()));
  return {HPX_MOVE(locality), first_index, last_index};
}

} // namespace hpx

namespace nw::graph {

template <>
struct graph_traits<usage_example::graph_type> {
  using vertex_id_type = default_vertex_id_type;
};

} // namespace nw::graph

#include "nwgraph/algorithms/partitioned_util.hpp"

namespace {

usage_example::graph_type make_partitioned_row_graph() {
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
  std::vector<vertex_id_t<std::remove_reference_t<Graph>>> result;
  auto row = static_cast<usage_example::neighbor_list>(G[u]);
  for (auto&& edge : row) {
    result.push_back(edge);
  }
  return result;
}

} // namespace

TEST_CASE("partitioned algorithms accept a partitioned_vector range-of-ranges graph",
          "[distributed][partitioned][concepts][smoke]") {
  static_assert(partition_token<usage_example::simple_partition>);
  static_assert(partitioned_graph<usage_example::graph_type>);
  static_assert(partitioned_algorithm_graph<usage_example::graph_type>);
  static_assert(std::same_as<partition_t<usage_example::graph_type>, usage_example::simple_partition>);

  auto G = make_partitioned_row_graph();

  // This example starts from hpx::partitioned_vector<std::vector<vertex>>.
  // The only distributed customization points it provides are partitions(G)
  // and vertex_partition(G, v); no custom graph local_view adapter is required.
  auto graph_partitions = partitions(G);
  REQUIRE(graph_partitions.size() == 2u);
  REQUIRE(vertex_partition(G, 0u) == graph_partitions[0]);
  REQUIRE(vertex_partition(G, 3u) == graph_partitions[1]);
  REQUIRE(is_local(graph_partitions[0], 0u));
  REQUIRE_FALSE(is_local(graph_partitions[0], 2u));
  REQUIRE(targets_of(G, 0u) == (std::vector<default_vertex_id_type>{1u, 2u}));
  REQUIRE(targets_of(G, 2u) == (std::vector<default_vertex_id_type>{0u, 1u, 3u}));

  auto avg_degree = partitioned_avg_degree_per_partition(G);
  REQUIRE(avg_degree == (std::vector<double>{2.0, 2.0}));

  auto avg_remote_degree = partitioned_avg_remote_degree_per_partition(G);
  REQUIRE(avg_remote_degree == (std::vector<double>{1.0, 1.0}));
}
