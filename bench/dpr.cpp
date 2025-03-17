/**
 * @file pr.cpp
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
 *
 */

static constexpr const char USAGE[] =
  R"(pr.exe: BGL17 page rank benchmark driver.
  Usage:
      pr.exe (-h | --help)
      pr.exe [--version ID...] -f FILE... [-i NUM] [-t NUM] [-n NUM] [-dvV] [--log FILE] [--log-header] [--partitions PARTS] [THREADS]...

  Options:
      -h, --help                show this screen
      --version ID              algorithm version to run [default: 0]
      -f FILE                   input file path
      -i NUM                    maximum iteration [default: 20]
      -t NUM                    tolerance [default: 1e-4]
      -n NUM                    number of trials [default: 1]
      --log FILE                log times to a file
      --log-header              add a header to the log file
      -d, --debug               run in debug mode
      -v, --verify              verify results
      -V, --verbose             run in verbose mode
      -p, --partitions PARTS    number of graph partitions to create [default: 1]
)";

#include <hpx/hpx_init.hpp>

#include <docopt.h>
#include "Log.hpp"
#include "common.hpp"

#include "nwgraph/partitioned_adjacency.hpp"
#include "nwgraph/algorithms/partitioned_page_rank_0.hpp"

#include <hpx/include/partitioned_vector.hpp>
#include <nwgraph/util/partitioned_serialize.hpp>

using unsigned_int = unsigned int;
HPX_REGISTER_PARTITIONED_VECTOR(unsigned_int)

HPX_REGISTER_PARTITIONED_VECTOR(float)


using namespace nw::graph::bench;
using namespace nw::graph;
using namespace nw::util;

template <typename Vector>
void print_n_ranks(const Vector& rankings, size_t n) {
  auto perm = proxysort<size_t>(rankings, std::greater<float>());
  for (size_t i = 0; i < 10; ++i) {
    std::cout << std::to_string(perm[i]) + ": " << std::to_string(rankings[perm[i]]) << std::endl;
  }
}

int hpx_main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  // Read the options
  bool verify = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug = args["--debug"].asBool();
  long trials = args["-n"].asLong() ? args["-n"].asLong() : 1;
  long max_iters = args["-i"].asLong() ? args["-i"].asLong() : 1;
  float tolerance = std::stof(args["-t"].asString());
  long num_partitions = args["--partitions"].asLong() ? args["--partitions"].asLong()
                                                      : hpx::get_num_localities(hpx::launch::sync);

  std::vector files = args["-f"].asStringList();
  std::vector ids = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  Times times;

  for (auto&& file : files) {
    //auto aos_a = load_graph<nw::graph::directedness::directed>(file);
    //if (verbose) {
    //  aos_a.stream_stats();
    //}

    //auto graph = build_adjacency<1>(aos_a);
    //if (verbose) {
    //  graph.stream_stats();
    //}

    //if (debug) {
    //  graph.stream_indices();
    //}

    if (!segment_files_exist(file, num_partitions)) {
      std::cout << "segment files do not exist, creating them\n";
      serialize_adjacency_graph_segments(file, num_partitions);
    }

    auto graph = load_partitioned_adjacency(file, num_partitions);

    auto degrees = build_degrees(graph);

    using vertex_id_type = typename decltype(graph)::vertex_id_type;
    using degree_type = typename decltype(degrees)::value_type;

    auto sizes = graph.indices_.get_partition_sizes();

    hpx::partitioned_vector<degree_type> p_degrees(
      degrees.size(), 0.0,
      hpx::explicit_container_layout(sizes, graph.indices_.get_partition_localities()));
    p_degrees.register_as("p_degrees");


    //hpx::copy(degrees.begin(), degrees.end(), p_degrees.begin());
    //degrees.register_as("p_degrees");
    for (size_t i = 0; i < graph.size(); ++i) {
      p_degrees[i] = degrees[i];
    }

    hpx::partitioned_vector<float> p_rankings(
      graph.indices_.size(), 0.0,
      hpx::explicit_container_layout(sizes, graph.indices_.get_partition_localities()));
    p_rankings.register_as("p_rankings");

    for (auto thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto id : ids) {
        for (size_t j = 0, e = trials; j < e; ++j) {
          times.record(
            file, id, thread,
            [&]
            {
              switch (id) {
              case 0:
                partitioned_page_rank_0(graph, p_degrees, p_rankings, 0.85f, tolerance, max_iters);
                break;

              default:
                std::cerr << "Unknown version id " << id << std::endl;
                break;
              }
            });
        }

        //if (verify) {
        //  std::cout << "Verifying\n";
        //  print_n_ranks(rankings, 10);
        //}
      }
    }
  }

  times.print(std::cout);

  if (args["--log"]) {
    auto file = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("pr", file, times, header, "Time(s)", "Tolerance");
  }

  hpx::finalize();
  return 0;
}

int main(int argc, char* argv[]) { return hpx::init(argc, argv); }
