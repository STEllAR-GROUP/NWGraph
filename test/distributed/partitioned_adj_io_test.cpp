
#include <iostream>
#include <queue>

#include "../common/test_header.hpp"


#include <hpx/hpx_init.hpp>
#include <hpx/include/partitioned_vector.hpp>

#include "nwgraph/graph_concepts.hpp"
#include "nwgraph/containers/aos.hpp"
#include "nwgraph/distributed/serialize.hpp"

// I don't like this being here
using unsigned_int = unsigned int;
HPX_REGISTER_PARTITIONED_VECTOR(unsigned_int)


using namespace nw::graph;
using namespace nw::util;

template <typename Graph>
bool contains(Graph& graph, size_t u, size_t v) {
  for (auto&& [x, y] : make_edge_range(graph)) {
    if (x == u && y == v) return true;
  }
  return false;
}

TEST_CASE("partitioned adj (compressed)  I/O", "[partitioned_compressed_io]") {
  SECTION("I/O (read real symmetric to edge_list and convert to compressed graph)") {
    //auto A = read_mm<directedness::directed>(DATA_DIR "tree.mmio");
    //auto B = read_mm<directedness::undirected>(DATA_DIR "tree.mmio");

    //auto C = read_mm<directedness::directed>(DATA_DIR "USAir97.mtx");
    //auto D = read_mm<directedness::undirected>(DATA_DIR "USAir97.mtx");

    auto A_local_edges = read_mm<directedness::directed>(DATA_DIR "karate.mtx");
    auto A_local = adjacency<0>(A_local_edges);

    
    std::string A_bin_file = partitioned_serialize_adj<0, directedness::directed>(DATA_DIR "karate.mtx");
    partitioned_adjacency A = partitioned_deserialize_adj<0, directedness::directed>(A_bin_file);

    REQUIRE(num_vertices(A) == num_vertices(A_local));
    REQUIRE(A.num_edges() == A_local.num_edges());

    for (auto [u, v] : make_edge_range(A)) {
      REQUIRE(contains(A_local, u, v));
    }

  }
  SECTION("I/O preserves adjacency<1> directed semantics") {
    REQUIRE_THROWS_AS(
      (partitioned_serialize_adj<1, directedness::directed>(DATA_DIR "karate.mtx")),
      std::logic_error);
    REQUIRE_THROWS_AS(
      (partitioned_deserialize_adj<1, directedness::directed>(DATA_DIR "karate.mtx")),
      std::logic_error);
  }
  SECTION("I/O preserves undirected adjacency semantics") {
    REQUIRE_THROWS_AS(
      (partitioned_serialize_adj<0, directedness::undirected>(DATA_DIR "karate.mtx")),
      std::logic_error);
    REQUIRE_THROWS_AS(
      (partitioned_deserialize_adj<0, directedness::undirected>(DATA_DIR "karate.mtx")),
      std::logic_error);
  }
  //SECTION("I/O (read pattern symmetric to edge_list and convert to compressed graph)") {
  //  auto A = read_mm<directedness::directed>(DATA_DIR "karate.mtx");
  //  auto B = read_mm<directedness::undirected>(DATA_DIR "karate.mtx");
  //}
  //SECTION("I/O (read real unsymmetric to edge_list and convert to compressed graph)") {
  //  auto A = read_mm<directedness::directed>(DATA_DIR "tree.mmio");
  //  auto B = read_mm<directedness::undirected>(DATA_DIR "tree.mmio");
  //}
  //SECTION("I/O (read pattern unsymmetric to edge_list and convert to compressed graph)") {}
  //SECTION("I/O (read to edge_list and convert to compressed matrix)") {}
}
