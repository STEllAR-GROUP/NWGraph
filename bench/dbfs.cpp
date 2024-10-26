/**
 * @file bfs.cpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * @authors
 *   Andrew Lumsdaine
 *   Tony Liu
 *   liux238
 *   Hartmut Kaiser
 *
 */

static constexpr char USAGE[] =
  R"(dbfs.exe: BGL17 breadth first search benchmark driver (distributed).
  Usage:
      dbfs.exe (-h | --help)
      dbfs.exe -f FILE [-r NODE | -s FILE] [-i NUM] [-a NUM] [-b NUM] [-B NUM] [-n NUM] [--seed NUM] [--version ID...] [--log FILE] [--log-header] [--partitions PARTS] [--batchsize SIZE] [-dvV] [THREADS]...

  Options:
      -h, --help              show this screen
      -f FILE                 input file path
      -i NUM                  number of iteration [default: 1]
      -a NUM                  alpha parameter [default: 15]
      -b NUM                  beta parameter [default: 18]
      -B NUM                  number of bins [default: 32]
      -n NUM                  number of trials [default: 1]
      -r NODE                 start from node r (default is random)
      -s, --sources FILE      sources file
      --seed NUM              random seed [default: 27491095]
      --version ID            algorithm version to run [default: 0]
      --log FILE              log times to a file
      --log-header            add a header to the log file
      -d, --debug             run in debug mode
      -v, --verify            verify results
      -V, --verbose           run in verbose mode
      -p, --partitions PARTS  number of graph partitions to create [default: 1]
      --batchsize SIZE        number asynchronous operations to batch [default: 10000]
)";

#ifndef NWGRAPH_HAVE_HPX
#error "This benchmark requires using the HPX backend for NWGraph"
#endif

#include <hpx/hpx_init.hpp>

#include "nwgraph/algorithms/bfs.hpp"
#include "nwgraph/algorithms/partitioned_bfs_1.hpp"
#include "nwgraph/partitioned_adjacency.hpp"

#include <filesystem>

#include <docopt.h>
#include "Log.hpp"
#include "common.hpp"

#include <hpx/include/partitioned_vector.hpp>

using unsigned_int = unsigned int;
HPX_REGISTER_PARTITIONED_VECTOR(unsigned_int)

using namespace nw::graph::bench;
using namespace nw::graph;
using namespace nw::util;

template <directedness Directedness, class... Attributes>
edge_list<Directedness, Attributes...> load_binary_graph(std::string file) {

  std::filesystem::path p(file), ext(".bmtk");
  p.replace_extension(ext);
  if (exists(p)) {
    edge_list<Directedness, Attributes...> el;
    el.deserialize(p.string());
    return el;
  }

  auto el = load_graph<Directedness, Attributes...>(file);
  el.serialize(p.string());
  return el;
}

template <adjacency_list_graph GraphT, class Vector>
auto edge_sizes(GraphT const& A, Vector const& vert_sizes) {

  Vector cedge_sizes;
  cedge_sizes.reserve(vert_sizes.size());

  auto begin = A.begin();
  auto prev_idx = begin.index();
  for (auto size : vert_sizes) {
    begin += size;
    cedge_sizes.push_back(begin.index() - prev_idx);
    prev_idx = begin.index();
  }

  return cedge_sizes;
}

auto partition_sizes(size_t num_partitions, size_t all_vertices) {

  std::vector<size_t> vert_sizes;
  vert_sizes.reserve(num_partitions);

  size_t part_size = (all_vertices + num_partitions - 1) / num_partitions;
  for (size_t part = 0, num_vertices = 0; part != num_partitions;
       ++part, num_vertices += part_size) {

    assert(all_vertices >= num_vertices);
    size_t this_part_size =
      (num_vertices + part_size > all_vertices ? all_vertices - num_vertices : part_size);

    vert_sizes.push_back(this_part_size);
  }

  return vert_sizes;
}

template <int Adj, directedness Directedness, typename... Attributes, class Vector>
auto build_partitioned_adjacency(edge_list<Directedness, Attributes...>& g,
                                 adjacency<Adj, Attributes...> const& A, Vector&& vert_sizes,
                                 Vector&& edge_sizes) {

  partitioned_adjacency<Adj, Attributes...> B(g, std::forward<Vector>(vert_sizes),
                                              std::forward<Vector>(edge_sizes), "local_pg",
                                              std::vector({hpx::find_here()}));

  std::copy(A.indices_.begin(), A.indices_.end(), B.indices_.begin());
  std::copy(A.to_be_indexed_.begin(), A.to_be_indexed_.end(), B.to_be_indexed_.begin());
  return B;
}

template <adjacency_list_graph adjacency_t>
void partitioned_copy(adjacency_t const& src, adjacency_t& dest) {

  dest.get_indices().copy_data_from(src.get_indices());
  dest.get_to_be_indexed().copy_data_from(src.get_to_be_indexed());
}


template <adjacency_list_graph Graph, class Vector>
auto distribute_compressed(size_t num_vertices, size_t num_edges, Graph& A, Vector&& vert_sizes,
                           Vector&& edge_sizes) {
  life_timer _(__func__);
  Graph B(num_vertices + 1, num_vertices + 1, num_edges, std::forward<Vector>(vert_sizes),
          std::forward<Vector>(edge_sizes), "pg", hpx::find_all_localities());
  partitioned_copy(A, B);
  return B;
}

int hpx_main(int argc, char* argv[]) {
  std::vector strings = std::vector<std::string>(argv + 1, argv + argc);
  std::map args = docopt::docopt(USAGE, strings, true);

  // Read the options
  bool verify = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug = args["--debug"].asBool();
  long trials = args["-n"].asLong() ? args["-n"].asLong() : 1;
  long iterations = args["-i"].asLong() ? args["-i"].asLong() : 1;
  long alpha = args["-a"].asLong() ? args["-a"].asLong() : 15;
  long beta = args["-b"].asLong() ? args["-b"].asLong() : 18;
  long num_bins = args["-B"].asLong() ? args["-B"].asLong() : 32;
  std::string file = args["-f"].asString();

  std::vector ids = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  long num_partitions = args["--partitions"].asLong() ? args["--partitions"].asLong()
                                                      : hpx::get_num_localities(hpx::launch::sync);
  long batchsize = args["--batchsize"].asLong() ? args["--batchsize"].asLong() : 10000;

  auto aos_a = load_binary_graph<nw::graph::directedness::directed>(file);

  if (verbose) {
    aos_a.stream_stats();
  }

  auto loc_graph = build_adjacency<1>(aos_a);
  auto gx = build_adjacency<0>(aos_a);

  auto cvert_sizes = partition_sizes(num_partitions, num_vertices(aos_a));
  auto cedge_sizes = edge_sizes(loc_graph, cvert_sizes);

  auto loc_part_graph = build_partitioned_adjacency(aos_a, loc_graph, cvert_sizes, cedge_sizes);

  auto graph = distribute_compressed(num_vertices(aos_a), num_edges(aos_a), loc_part_graph,
                                     std::move(cvert_sizes), std::move(cedge_sizes));

  // free non-needed memory
  loc_graph = adjacency<1>{};
  loc_part_graph = partitioned_adjacency<1>{};

  if (verbose) {
    graph.stream_stats();
  }

  if (debug) {
    graph.stream_indices();
  }

  using vertex_id_type = vertex_id_t<decltype(graph)>;

  std::vector<vertex_id_type> sources;
  if (args["--sources"]) {
    sources = load_sources_from_file(graph, args["--sources"].asString());
    trials = sources.size();
  }
  else if (args["-r"]) {
    sources.resize(trials);
    std::fill(sources.begin(), sources.end(), args["-r"].asLong());
  }
  else {
    sources = build_random_sources(graph, trials, args["--seed"].asLong());
  }

  Times<vertex_id_type> times;

  std::map<long, std::vector<size_t>> levels;

  for (auto&& thread : threads) {
    auto _ = set_n_threads(thread);
    for (auto&& id : ids) {
      for (auto&& source : sources) {
        if (verbose) {
          std::cout << "source: " << source << "\n";
        }

        auto&& [time, parents] = time_op(
          [&]
          {
            switch (id) {
            case 0:
              return bfs(graph, source);

            case 1:
              return partitioned_bfs_1(graph, source, batchsize);
#if 0
            case 1:
              return bfs_v1(graph, gx, source, num_bins, alpha, beta);
            case 2:
              return bfs_v2(graph, gx, source, num_bins, alpha, beta);
            case 6:
              return bfs_v6(graph, source);
            case 7:
              return bfs_v7(graph, source);
            // case 8:
            //   return bfs_v8(graph, source);
            // case 9:
            //   return bfs_v9(graph, source);
            // case 10:
            //   return bfs_top_down(graph, source);
            case 11:
              return bfs(graph, gx, source, num_bins, alpha, beta);
            // case 12:
            //   return bfs_top_down_bitmap(graph, source);
            // case 13:
            //   return bfs_bottom_up(graph, gx, source);
#endif
            default:
              std::cerr << "Unknown version " << id << "\n";
              return std::vector<vertex_id_type>();
            }
          });

        if (verify) {
          BFSVerifier(graph, gx, source, parents);
        }

        times.append(file, id, thread, time, source);
      }
    }
  }

  times.print(std::cout);

  if (args["--log"]) {
    auto file = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("bfs", file, times, header, "Time(s)", "Source");
  }

  hpx::finalize();
  return 0;
}

int main(int argc, char* argv[]) { return hpx::init(argc, argv); }
