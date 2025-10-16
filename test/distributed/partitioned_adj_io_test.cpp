
#include <iostream>
#include <queue>

#include "../common/test_header.hpp"


#include <hpx/hpx_init.hpp>
#include <hpx/include/partitioned_vector.hpp>

#include "nwgraph/containers/aos.hpp"
#include "nwgraph/util/partitioned_serialize.hpp"

// I don't like this being here
using unsigned_int = unsigned int;
HPX_REGISTER_PARTITIONED_VECTOR(unsigned_int)


using namespace nw::graph;
using namespace nw::util;

template <typename Graph>
bool contains(Graph graph, size_t u, size_t v) {
  for (auto&& [x, y] : graph) {
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

    auto A_local = read_mm<directedness::undirected>(DATA_DIR "karate.mtx");

    
    std::string A_bin_file = partitioned_serialize_adj(DATA_DIR "karate.mtx");
    auto A = partitioned_deserialize_adj(A_bin_file);

    REQUIRE(A.num_vertices() == A_local.num_vertices());
    REQUIRE(A.num_edges() == A_local.num_edges());


    //for (auto rng : A) {
    //    for (auto edge : rng) {
    //        auto u = std::get<0>(edge);
    //      auto v = std::get<1>(edge);
    //    REQUIRE(contains(A_local, u, v));
    //    }
    //}

    for (auto [u, v] : make_edge_range(A)) {
      REQUIRE(contains(A_local, u, v));
    }

   



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
