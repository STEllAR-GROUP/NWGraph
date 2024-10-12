/**
 * @file dtc.cpp
 *
 * @copyright Copyright 2014, Software Engineering Institute, Carnegie Mellon University
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

static constexpr char USAGE[] =
  R"(dtc.exe: BGL17 triangle counting benchmark driver (distributed).
  Usage:
      dtc.exe (-h | --help)
      dtc.exe -f FILE... [--version ID...] [-n NUM] [--lower | --upper] [--relabel] [--heuristic] [--log FILE] [--log-header] [--format FORMAT] [--partitions PARTS] [-dvV] [THREADS...]

  Options:
      -h, --help              show this screen
      --version ID            algorithm version to run [default: 4]
      -f FILE                 input file path
      -n NUM                  number of trials [default: 1]
      --lower                 lower triangular order [default: false]
      --upper                 upper triangular order [default: true]
      --relabel               relabel the graph
      --heuristic             use heuristic to decide whether to relabel or not
      --format FORMAT         specify which graph storage format [default: CSR]
      --log FILE              log times to a file
      --log-header            add a header to the log file
      -d, --debug             run in debug mode
      -v, --verify            verify results
      -V, --verbose           run in verbose mode
      -p, --partitions PARTS  number of graph partitions to create [default: 1]
)";


#ifndef NWGRAPH_HAVE_HPX
#error "This benchmark requires using the HPX backend for NWGraph"
#endif

#include <hpx/hpx_init.hpp>

#include "nwgraph/adjacency.hpp"
#include "nwgraph/edge_list.hpp"
#include "nwgraph/volos.hpp"
#include "nwgraph/vovos.hpp"

#include <docopt.h>
#include <filesystem>
#include <tuple>
#include "common.hpp"
#include "nwgraph/algorithms/partitioned_triangle_count_0.hpp"
#include "nwgraph/algorithms/partitioned_triangle_count_1.hpp"
#include "nwgraph/algorithms/partitioned_triangle_count_2.hpp"
#include "nwgraph/partitioned_adjacency.hpp"

#include <date/date.h>
#include "config.h"

#include <nlohmann/json.hpp>

#include <hpx/include/partitioned_vector.hpp>

using unsigned_int = unsigned int;
HPX_REGISTER_PARTITIONED_VECTOR(unsigned_int)

using json = nlohmann::json;

using namespace nw::graph::bench;
using namespace nw::graph;
using namespace nw::util;

template <class Vector>
static void tc_relabel(edge_list<directedness::undirected>& A, Vector&& degrees,
                       std::string const& direction) {
  life_timer _(__func__);
  relabel_by_degree<0>(A, direction, degrees);
}

template <size_t id = 0>
static void clean(edge_list<directedness::undirected>& A, std::string const& succession) {
  life_timer _(__func__);
  swap_to_triangular<id>(A, succession);
  lexical_sort_by<id>(A);
  uniq(A);
  remove_self_loops(A);
}

auto vertex_sizes(size_t num_partitions, size_t all_vertices) {

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

inline auto compress(edge_list<directedness::undirected>& A) {
  life_timer _(__func__);
  adjacency<0> B(num_vertices(A));
  push_back_fill(A, B);
  return B;
}

template <adjacency_list_graph adjacency_t, typename... Ts>
auto partitioned_push_back_fill_helper(size_t idx, adjacency_t& cs,
                                       std::tuple<Ts...> const& theTuple) {
  std::apply([&](Ts const&... args) { cs.push_at(idx, args...); }, theTuple);
}

template <edge_list_c edge_list_t, adjacency_list_graph adjacency_t>
void partitioned_push_back_fill(edge_list_t& el, adjacency_t& cs) {
  cs.open_for_push_back();

  std::for_each(el.begin(), el.end(), [&, idx = 0](auto&& elt) mutable
                { partitioned_push_back_fill_helper(idx++, cs, elt); });

  cs.close_for_push_back();
}

template <adjacency_list_graph adjacency_t>
void partitioned_copy(adjacency_t const& src, adjacency_t& dest) {

  dest.get_indices().copy_data_from(src.get_indices());
  dest.get_to_be_indexed().copy_data_from(src.get_to_be_indexed());
}

template <adjacency_list_graph Graph, class Vector>
auto compress(edge_list<directedness::undirected>& A, Vector&& vert_sizes, Vector&& edge_sizes) {
  life_timer _(__func__);
  Graph B(num_vertices(A), num_edges(A), std::forward<Vector>(vert_sizes),
          std::forward<Vector>(edge_sizes), "local_pg", std::vector({hpx::find_here()}));
  partitioned_push_back_fill(A, B);
  return B;
}

template <adjacency_list_graph Graph, class Vector>
auto distribute_compressed(size_t num_vertices, size_t num_edges, Graph& A, Vector&& vert_sizes,
                           Vector&& edge_sizes) {
  life_timer _(__func__);
  Graph B(num_vertices, num_edges, std::forward<Vector>(vert_sizes),
          std::forward<Vector>(edge_sizes), "pg", hpx::find_all_localities());
  partitioned_copy(A, B);
  return B;
}

// heuristic to see if sufficiently dense power-law graph
template <edge_list_graph EdgeList, class Vector>
static bool worth_relabeling(EdgeList const& el, Vector const& degree) {
  using vertex_id_type = typename EdgeList::vertex_id_type;

  int64_t average_degree = el.size() / (el.num_vertices()[0]);
  if (average_degree < 10)
    return false;

  int64_t num_samples = std::min<int64_t>(1000L, el.num_vertices()[0]);
  int64_t sample_total = 0;
  std::vector<int64_t> samples(num_samples);

  std::mt19937 rng;
  std::uniform_int_distribution<vertex_id_type> udist(0, el.num_vertices()[0] - 1);

  for (int64_t trial = 0; trial < num_samples; trial++) {
    samples[trial] = degree[udist(rng)];
    sample_total += samples[trial];
  }
#ifdef NWGRAPH_HAVE_HPX
  hpx::sort(hpx::execution::par_unseq, samples.begin(), samples.end());
#else
  std::sort(std::execution::par_unseq, samples.begin(), samples.end());
#endif
  double sample_average = static_cast<double>(sample_total) / num_samples;
  double sample_median = samples[num_samples / 2];
  return sample_average / 1.3 > sample_median;
}

// Taken from GAP and adapted to NW Graph
template <adjacency_list_graph Graph>
static size_t TCVerifier(Graph& graph) {
  using vertex_id_type = typename Graph::vertex_id_type;

  life_timer _(__func__);
  std::vector<std::tuple<vertex_id_type>> intersection;
  intersection.reserve(graph.size());
  for (auto&& [u, v] : edge_range(graph)) {
    auto u_out = graph[u];
    auto v_out = graph[v];
    std::set_intersection(u_out.begin(), u_out.end(), v_out.begin(), v_out.end(),
                          std::back_inserter(intersection));
  }
  size_t total = intersection.size();
  return total; // note that our processed Graph doesn't produce extra counts
                // like the GAP verifier normally would
}

auto config_log() {
  std::string uuid_;
  char host_[16];
  std::string date_;
  std::string git_branch_;
  std::string git_version_;
  size_t uuid_size_ = 24;

  auto seed = std::random_device();
  auto gen = std::mt19937(seed());
  auto dis = std::uniform_int_distribution<short>(97, 122);
  uuid_.resize(uuid_size_);
  std::generate(uuid_.begin(), uuid_.end(), [&] { return dis(gen); });

  if (int e = gethostname(host_, sizeof(host_))) {
    std::cerr << "truncated host name\n";
    strncpy(host_, "ghost", 15);
  }
  {
    std::stringstream ss;
    date::year_month_day date = date::floor<date::days>(std::chrono::system_clock::now());
    ss << date;
    date_ = ss.str();
  }

  if (!std::system("git rev-parse --abbrev-ref HEAD > git_branch.txt")) {
    std::ifstream("git_branch.txt") >> git_branch_;
  }
  if (!std::system("git log --pretty=format:'%h' -n 1 > git_version.txt")) {
    std::ifstream("git_version.txt") >> git_version_;
    if (std::system("git diff --quiet --exit-code")) {
      git_version_ += "+";
    }
  }

  json config = {{"Host", host_},
                 {"Date", date_},
                 {"git_branch", git_branch_},
                 {"git_version", git_version_},
                 {"Build", BUILD_TYPE},
                 {"CXX_COMPILER", CXX_COMPILER},
                 {"CXX_COMPILER_ID", CXX_COMPILER_ID},
                 {"CXX_VERSION", CXX_VERSION}};

  return config;
}

template <typename Args>
auto args_log(Args const& args) {
  json arg_log;

  for (auto&& arg : args) {
    std::stringstream buf;
    buf << std::get<1>(arg);
    arg_log.push_back({std::get<0>(arg), buf.str()});
  }
  return arg_log;
}

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

template <typename Graph>
void run_bench(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  // Read the easy options
  bool verify = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug = args["--debug"].asBool();
  long trials = args["-n"].asLong() ? args["-n"].asLong() : 1;
  long num_partitions = args["--partitions"].asLong() ? args["--partitions"].asLong()
                                                      : hpx::get_num_localities(hpx::launch::sync);

  // Read the more complex options
  std::string direction = "ascending";
  std::string succession = "successor";
  if (args["--lower"].asBool()) {
    direction = "descending";
    succession = "predecessor";
  }

  std::vector files = args["-f"].asStringList();
  std::vector ids = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  json file_log = {};
  size_t file_ctr = 0;
  for (auto&& file : files) {
    std::cout << "processing " << file << "\n";

    auto el_a = load_binary_graph<nw::graph::directedness::undirected>(file);
    auto degree = degrees(el_a);

    // Run and time relabeling. This operates directly on the incoming edgelist.
    bool relabeled = false;
    auto&& [relabel_time] = time_op(
      [&]
      {
        if (args["--relabel"].asBool()) {
          if (args["--heuristic"].asBool() == false || worth_relabeling(el_a, degree)) {
            tc_relabel(el_a, degree, direction);
            relabeled = true;
          }
        }
      });

    // Force relabel time to 0 if we didn't relabel, this just suppresses
    // nanosecond noise from the time_op function when relabeling wasn't done.
    if (!relabeled) {
      relabel_time = 0.0;
    }

    // Clean up the edgelist to deal with the normal issues related to
    // undirectedness.
    auto&& [clean_time] = time_op([&] { clean<0>(el_a, succession); });

    auto local_cel_a = compress(el_a);

    auto cvert_sizes = vertex_sizes(num_partitions, num_vertices(el_a));
    auto cedge_sizes = edge_sizes(local_cel_a, cvert_sizes);

    // create a local copy of the partitioned graph
    auto l_cel_a = compress<Graph>(el_a, cvert_sizes, cedge_sizes);

    // now distribute the partitioned data
    auto cel_a = distribute_compressed(num_vertices(el_a), num_edges(el_a), l_cel_a,
                                       std::move(cvert_sizes), std::move(cedge_sizes));

    //    if (debug) {
    // cel_a.stream_indices();
    //}

    // If we're verifying then compute the number of triangles once for this
    // graph.
    size_t v_triangles = 0;
    if (verify) {
      v_triangles = TCVerifier(cel_a);
      std::cout << "verifier reports " << v_triangles << " triangles\n";
    }

    json thread_log = {};
    size_t thread_ctr = 0;

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);

      json id_log = {};
      size_t id_ctr = 0;
      for (auto&& id : ids) {

        json run_log = {};
        size_t run_ctr = 0;

        for (int j = 0; j < trials; ++j) {
          if (verbose) {
            std::cout << "running version:" << id << " threads:" << thread << "\n";
          }

          auto&& [time, triangles] = time_op(
            [&]() -> size_t
            {
              switch (id) {
              case 0:
                return triangle_count(cel_a);

              case 1:
                return partitioned_triangle_count_0(cel_a);

              case 2:
                return partitioned_triangle_count_1(cel_a);

              case 3:
                return partitioned_triangle_count_2(cel_a);
#if 0
              case 1:
                return triangle_count_v1(cel_a);
              case 2:
                return triangle_count_v2(cel_a);
              case 3:
                return triangle_count_v3(cel_a);
              case 4:
                return triangle_count(cel_a, thread);
               case 5:
                 return triangle_count_v5(cel_a.begin(), cel_a.end(), thread);
               case 6:
                 return triangle_count_v6(cel_a.begin(), cel_a.end(), thread);
               case 7:
                 return triangle_count_v7(cel_a);
               case 8:
                 return triangle_count_v7(cel_a, std::execution::seq, std::execution::par_unseq);
               case 9:
                 return triangle_count_v7(cel_a, std::execution::par_unseq, std::execution::par_unseq);
               case 10:
                 return triangle_count_v10(cel_a);
               case 11:
                 return triangle_count_v10(cel_a, std::execution::par_unseq, std::execution::par_unseq, std::execution::par_unseq);
               case 12:
                 return triangle_count_v12(cel_a, thread);
               case 13:
                 return triangle_count_v13(cel_a, thread);
               case 14:
                 return triangle_count_v14(cel_a);
               case 15:
                 return triangle_count_edgesplit(cel_a, thread);
               case 16:
                 return triangle_count_edgesplit_upper(cel_a, thread);
#ifdef ONE_DIMENSIONAL_EDGE
               case 17:
                 return triangle_count_edgerange(cel_a);
               case 18:
                 return triangle_count_edgerange_cyclic(cel_a, thread);
#endif
#endif
              default:
                std::cerr << "Unknown version id " << id << "\n";
                return 0ul;
              }
            });

          run_log[run_ctr++] = {{"id", id},
                                {"num_threads", thread},
                                {"trial", j},
                                {"elapsed", time},
                                {"elapsed+relabel", time + relabel_time},
                                {"triangles", triangles}};

          if (verify && triangles != v_triangles) {
            std::cerr << "Inconsistent results: v" << id << " failed verification for " << file
                      << " using " << thread << " threads (reported " << triangles << ")\n";
          }
        } // for j in trials

        id_log[id_ctr++] = {{"id", id}, {"runs", std::move(run_log)}};
      } // for id in ids

      thread_log[thread_ctr++] = {{"num_thread", thread}, {"runs", std::move(id_log)}};
    } // for thread in threads

    file_log[file_ctr++] = {
      {"File", file},           {"Relabel_time", relabel_time}, {"Clean_time", clean_time},
      {"Relabeled", relabeled}, {"Num_trials", trials},         {"Runs", std::move(thread_log)}};

  } // for each file
  if (args["--log"]) {

    json log_log = {
      {"Config", config_log()}, {"Args", args_log(args)}, {"Files", std::move(file_log)}};

    if (args["--log"].asString() == "-") {
      std::cout << log_log << std::endl;
    }
    else {
      std::ofstream outfile(args["--log"].asString(), std::ios_base::app);
      outfile << log_log << std::endl;
    }
  }
}

int hpx_main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  if (args["--format"].asString() == "CSR") {
    run_bench<partitioned_adjacency<0>>(argc, argv);

    //  } else if (args["--format"].asString() == "VOV") {
    //    run_bench<vov<0>>(argc, argv);
    //
    //  } else if (args["--format"].asString() == "VOL") {
    //    run_bench<adj_list<0>>(argc, argv);
  }
  else {
    std::cerr << "bad format" << std::endl;
    hpx::finalize();
    return -1;
  }

  hpx::finalize();
  return 0;
}

int main(int argc, char* argv[]) { return hpx::init(argc, argv); }
