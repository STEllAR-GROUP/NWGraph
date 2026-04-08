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
#include "nwgraph/distributed/algorithms/bfs_1.hpp"
#include "nwgraph/distributed/adjacency.hpp"
#include "nwgraph/distributed/algorithms/util.hpp"
#include <nwgraph/distributed/serialize.hpp>

#include <hpx/include/partitioned_vector.hpp>

#include <filesystem>

#include <docopt.h>
#include "Log.hpp"
#include "common.hpp"


using unsigned_int = unsigned int;
HPX_REGISTER_PARTITIONED_VECTOR(unsigned_int)

using namespace nw::graph::bench;
using namespace nw::graph;
using namespace nw::util;

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

  auto graph = load_partitioned_adjacency_graph<0, directedness::directed>(file);

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
            default:
              std::cerr << "Unsupported distributed BFS version id " << id
                        << "; available versions: 0, 1\n";
              return std::vector<vertex_id_type>();
            }
          });

        if (verify) {
          auto aos_a = load_binary_graph<nw::graph::directedness::directed>(file);
          auto loc_graph = build_adjacency<1>(aos_a);
          auto gx = build_adjacency<0>(aos_a);
          BFSVerifier(loc_graph, gx, source, parents);
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
