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
      pr.exe [--version ID...] -f FILE... [-i NUM] [-t NUM] [-n NUM] [-dvV] [--log FILE] [--log-header] [--partitions PARTS] [--batchsize SIZE] [THREADS]...

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
      -p, --partitions PARTS    number of graph partitions to create [default: 0]
      -b, --batchsize SIZE    number asynchronous operations to batch [default: 10000]
)";

#include <hpx/hpx_init.hpp>

#include <docopt.h>
#include "Log.hpp"
#include "common.hpp"

#include <nwgraph/util/partitioned_serialize.hpp>
#include "nwgraph/algorithms/partitioned_page_rank_0.hpp"
#include "nwgraph/algorithms/partitioned_page_rank_1.hpp"
#include "nwgraph/algorithms/partitioned_page_rank_2.hpp"
#include "nwgraph/algorithms/partitioned_page_rank_3.hpp"
#include "nwgraph/algorithms/partitioned_page_rank_4.hpp"
#include "nwgraph/algorithms/partitioned_util.hpp"
#include "nwgraph/experimental/algorithms/page_rank.hpp"
#include "nwgraph/partitioned_adjacency.hpp"


#include <hpx/include/partitioned_vector.hpp>
#include <nwgraph/partitioned_build.hpp>

using unsigned_int = unsigned int;
HPX_REGISTER_PARTITIONED_VECTOR(unsigned_int)

HPX_REGISTER_PARTITIONED_VECTOR(float)


using namespace nw::graph::bench;
using namespace nw::graph;
using namespace nw::util;

namespace {

constexpr float kVerificationRelTolerance = 0.001f;
constexpr float kVerificationAbsTolerance = 1.0e-6f;

template <typename T>
std::vector<T> copy_from_partitioned_vector(hpx::partitioned_vector<T>& src) {
  auto const sizes = src.get_partition_sizes();
  auto const src_partitions = sizes.size();
  std::vector<std::size_t> const empty;
  std::vector<T> local_values;
  local_values.reserve(src.size());
  for (std::size_t part = 0; part != src_partitions; ++part) {
    auto values = src.get_values(hpx::launch::sync, part, empty);
    local_values.insert(local_values.end(), values.begin(), values.end());
  }
  return local_values;
}

bool verify_partitioned_pagerank(
  const std::string& file, hpx::partitioned_vector<float>& distributed_rankings, float tolerance,
  long max_iters) {
  auto edge_list = load_binary_graph<nw::graph::directedness::directed>(file);
  auto local_graph = build_adjacency<1>(edge_list);
  auto local_degrees = build_degrees(local_graph);
  std::vector<float> local_rankings(local_graph.size());

  {
    nw::util::life_timer _("local verification pagerank");
    page_rank_v1(local_graph, local_degrees, local_rankings, 0.85f, tolerance, max_iters);
  }

  auto gathered_rankings = copy_from_partitioned_vector(distributed_rankings);
  float max_abs_error = 0.0f;
  float max_rel_error = 0.0f;

  for (std::size_t i = 0; i < gathered_rankings.size(); ++i) {
    float distributed = gathered_rankings[i];
    float reference = local_rankings[i];
    float greater = std::max(std::fabs(distributed), std::fabs(reference));
    float abs_error = std::fabs(distributed - reference);
    float rel_error = greater > 0.0f ? abs_error / greater : abs_error;
    float compare_tolerance = std::max(kVerificationAbsTolerance, kVerificationRelTolerance * greater);

    max_abs_error = std::max(max_abs_error, abs_error);
    max_rel_error = std::max(max_rel_error, rel_error);

    if (abs_error >= compare_tolerance) {
      std::cerr << "Results do not match\n";
      std::cerr << "First mismatch at vertex " << i << ": distributed=" << distributed
                << ", reference=" << reference << ", abs_error=" << abs_error
                << ", max_abs_error=" << max_abs_error
                << ", max_rel_error=" << max_rel_error << "\n";
      return false;
    }
  }

  std::cerr << "Verification passed with max_abs_error=" << max_abs_error
            << " and max_rel_error=" << max_rel_error << "\n";
  return true;
}

} // namespace

int hpx_main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  // Read the options
  bool verify = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug = args["--debug"].asBool();
  long trials = args["-n"].asLong();
  long max_iters = args["-i"].asLong();
  float tolerance = std::stof(args["-t"].asString());
  long num_partitions = args["--partitions"].asLong() ? args["--partitions"].asLong()
                                                      : hpx::get_num_localities(hpx::launch::sync);
  long batchsize = args["--batchsize"].asLong();

  std::cout << "Running on " << num_partitions << " partitions\n";

  std::vector files = args["-f"].asStringList();
  std::vector ids = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  Times<float> times;

  for (auto&& file : files) {

    partitioned_adjacency graph = load_partitioned_adjacency_graph(file);

    auto sizes = graph.indices_.get_partition_sizes();

    auto p_degrees = partitioned_degrees(graph);

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
              case 1:
                partitioned_page_rank_1(graph, p_degrees, p_rankings, 0.85f, tolerance, max_iters,
                                        batchsize);
                break;
              case 2:
                // This doesn't need a seperate degrees vector, because it assumes the inverse
                // edge directionality, meaning that the out-degrees is the size of the adjacency
                // list of each vertex.
                partitioned_page_rank_2(graph, p_rankings, 0.85f, tolerance, max_iters);
                break;
              case 3:
                partitioned_page_rank_3(graph, p_degrees, p_rankings, 0.85f, tolerance, max_iters,
                                        batchsize);
                break;
              case 4: 
                  partitioned_page_rank_4(graph, p_degrees, p_rankings, 0.85f, tolerance, max_iters,
                                        batchsize);
                break;
              default:
                std::cerr << "Unknown version id " << id << std::endl;
                break;
              }
            },
            tolerance);
        }

        if (verify && id == 2) {
          std::cout << "Skipping verification for partitioned_page_rank_2, as it assumes inverse "
                       "edge directionality."
                    << std::endl;
        }
        else if (verify) {
          std::cout << "Verifying..." << std::endl;
          nw::util::life_timer _("verification");
          verify_partitioned_pagerank(file, p_rankings, tolerance, max_iters);
        }
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
