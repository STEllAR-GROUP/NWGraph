#include <filesystem>
#include <fstream>
#include <iostream>

#include "nwgraph/adjacency.hpp"
#include "nwgraph/algorithms/partitioned_algorithm.hpp"
#include "nwgraph/edge_list.hpp"
#include "nwgraph/graph_base.hpp"
#include "nwgraph/io/mmio.hpp"
#include "nwgraph/partitioned_adjacency.hpp"

#include <hpx/async_combinators/wait_all.hpp>
#include <hpx/executors/execution_policy.hpp>
#include <hpx/include/partitioned_vector_predef.hpp>
#include <hpx/parallel/segmented_algorithms/detail/dispatch.hpp>


namespace nw::graph {

  namespace detail {

    auto partitioned_vertex_sizes(size_t num_partitions, size_t all_vertices) {

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

    auto compress_partial(edge_list<directedness::undirected>& A, size_t n_indices,
                          size_t n_to_be_indexed, size_t idx_offset) {

      std::vector<default_index_t> indices(n_indices, 0);
      std::vector<default_index_t> to_be_indexed(n_to_be_indexed);

      size_t idx_start = std::get<0>(A[0]);
      for (size_t i = 0; i < A.size(); i++) {
        indices[std::get<0>(A[i]) - idx_start]++;
        to_be_indexed[i] = std::get<1>(A[i]);
      }

      std::exclusive_scan(indices.begin(), indices.end(), indices.begin(), idx_offset);

      return std::make_tuple(std::move(indices), std::move(to_be_indexed));
    };

    std::string get_adj_filename(std::string mtx_file) {
      return mtx_file.substr(0, mtx_file.find_last_of('.')) + ".adj";
    }

  } // namespace detail

  // max_part_size limits the size of the .mtx file that is loaded into memory at once,
  // which is useful for very large graphs.
  void serialize_adj(std::string mtx_file, size_t max_part_size = 2 << 26) {

    std::string file_name = detail::get_adj_filename(mtx_file);

    if (std::filesystem::exists(file_name)) {
      std::cout << "Adjacency file already exists: " << file_name << ". Skipping serialization."
                << std::endl;
      return;
    }

    // Get matrix size
    std::ifstream in_stream(mtx_file);
    auto [n_vertices, _, n_edges] = read_mm_metadata(in_stream);
    size_t n_parts = (n_vertices + max_part_size - 1) / max_part_size;
    auto part_sizes = detail::partitioned_vertex_sizes(n_parts, n_vertices);


    char magic[] = "NWGRAPH ADJACENCY BINARY FILE";

    size_t header_size = sizeof(magic) + sizeof(n_vertices) + sizeof(n_edges);
    size_t indices_size = sizeof(unsigned int) * (n_vertices + 1); // for adjacency list
    size_t to_be_indexed_size = sizeof(unsigned int) * (n_edges); // for to_be_indexed
    size_t total_size = header_size + indices_size + to_be_indexed_size;
    // Create file
    {
      std::ofstream f(file_name /*, std::ios::binary*/);
    }
    std::filesystem::resize_file(file_name,
                                 total_size); // Resize file to accommodate header and data

    // Write header
    {
      std::ofstream f(file_name, std::ios::binary);
      f.write(magic, sizeof(magic));
      f.write(reinterpret_cast<const char*>(&n_vertices), sizeof(n_vertices));
      f.write(reinterpret_cast<const char*>(&n_edges), sizeof(n_edges));
      assert(f.tellp() == header_size);
    }


    // Open file twice for writing indices and to_be_indexed simultaneously
    auto open_mode = std::ios::binary | std::ios::in | std::ios::out | std::ios::ate;
    std::ofstream f_idx(file_name, open_mode);
    f_idx.seekp(header_size); // Move to the end of the header
    std::ofstream f_to_be_idx(file_name, open_mode);
    f_to_be_idx.seekp(header_size + indices_size); // Move to the end of the indices

    // Keep count for sanity
    size_t total_written_indices = 0;
    size_t total_written_to_be_indexed = 0;
    size_t idx_start = 0;
    for (size_t part_size : part_sizes) {
      size_t idx_end = idx_start + part_size;

      // Reads the whole file, but only keeps edges whose source index is within the range
      auto pred = [idx_start, idx_end](auto&& d0, auto&& d1, auto&& v)
      { return d0 >= idx_start && d0 < idx_end; };

      in_stream.seekg(0); // Reset to beggining of the file
      auto edgelist = read_mm<decltype(pred), directedness::undirected>(in_stream, pred);

      // Create partial adjacency
      sort_by<0>(edgelist);
      auto [indices, to_be_indexed] = detail::compress_partial(
        edgelist, idx_end - idx_start, edgelist.size(), total_written_to_be_indexed);

      auto f_idx_pos = f_idx.tellp();

      // Serialize the indices and to_be_indexed vectors
      f_idx.write(reinterpret_cast<const char*>(indices.data()),
                  indices.size() * sizeof(unsigned int));
      f_idx.flush();

      assert((f_idx.tellp() - f_idx_pos) / sizeof(unsigned int) == indices.size());

      assert(f_idx.tellp() <= header_size + indices_size);

      auto f_to_be_idx_pos = f_to_be_idx.tellp();

      f_to_be_idx.write(reinterpret_cast<const char*>(to_be_indexed.data()),
                        to_be_indexed.size() * sizeof(unsigned int));
      f_to_be_idx.flush();

      assert((f_to_be_idx.tellp() - f_to_be_idx_pos) / sizeof(unsigned int) ==
             to_be_indexed.size());

      assert(f_to_be_idx.tellp() <= total_size);

      total_written_indices += indices.size();
      total_written_to_be_indexed += to_be_indexed.size();
      idx_start = idx_end; // Update start index for next partition
    }
    assert(total_written_indices == n_vertices);
    assert(total_written_to_be_indexed == n_edges);
    // Write final idx, which should point to the end of the last to_be_indexed
    f_idx.write(reinterpret_cast<const char*>(&n_edges), sizeof(unsigned int));
    // final checks
    assert(f_idx.tellp() == header_size + indices_size);
    assert(f_to_be_idx.tellp() == total_size);

    f_idx.close();
    f_to_be_idx.close();
    std::cout << "Adjacency file serialized successfully." << std::endl;
    std::cout << "File size: " << total_size << " bytes." << std::endl;
  }


  auto deserialize_adj_part(std::string mtx_file, size_t begin_idx, size_t end_idx) {
    assert(begin_idx <= end_idx);
    std::string file_name = detail::get_adj_filename(mtx_file);
    std::ifstream f(file_name, std::ios::binary | std::ios::in);
    char magic[30];
    f.read(magic, 30);
    if (strncmp(magic, "NWGRAPH ADJACENCY BINARY FILE", 30) != 0) {
      std::cerr << "Error: Invalid file format\n";
      exit(1);
    }
    size_t n_vertices;
    size_t n_edges;

    f.read(reinterpret_cast<char*>(&n_vertices), sizeof(n_vertices));
    f.read(reinterpret_cast<char*>(&n_edges), sizeof(n_edges));

    assert(end_idx <= n_vertices);

    size_t header_size = sizeof(magic) + sizeof(n_vertices) + sizeof(n_edges);
    size_t indices_size = sizeof(unsigned int) * (n_vertices + 1); // for adjacency list
    size_t to_be_indexed_size = sizeof(unsigned int) * (n_edges); // for to_be_indexed

    std::vector<unsigned int> indices(end_idx - begin_idx + 1);
    f.seekg(header_size + begin_idx * sizeof(unsigned int));
    f.read(reinterpret_cast<char*>(indices.data()), indices.size() * sizeof(unsigned int));

    size_t begin_to_be_idx = indices.front();
    size_t end_to_be_idx = indices.back();

    std::vector<unsigned int> to_be_indexed(end_to_be_idx - begin_to_be_idx);
    f.seekg(header_size + indices_size + begin_to_be_idx * sizeof(unsigned int));
    f.read(reinterpret_cast<char*>(to_be_indexed.data()),
           to_be_indexed.size() * sizeof(unsigned int));

    return std::make_tuple(std::move(indices), std::move(to_be_indexed));
  }

  std::vector<size_t> read_edge_sizes(std::string file_name, std::vector<size_t> idx_sizes) {
    std::ifstream f(file_name, std::ios::binary | std::ios::in);
    char magic[30];
    f.read(magic, 30);
    if (strncmp(magic, "NWGRAPH ADJACENCY BINARY FILE", 30) != 0) {
      std::cerr << "Error: Invalid file format\n";
      exit(1);
    }
    size_t n_vertices;
    size_t n_edges;
    f.read(reinterpret_cast<char*>(&n_vertices), sizeof(n_vertices));
    f.read(reinterpret_cast<char*>(&n_edges), sizeof(n_edges));
    size_t header_size = sizeof(magic) + sizeof(n_vertices) + sizeof(n_edges);
    size_t indices_size = sizeof(unsigned int) * (n_vertices + 1); // for adjacency list
    size_t to_be_indexed_size = sizeof(unsigned int) * (n_edges); // for to_be_indexed
    f.seekg(header_size + indices_size);
    std::vector<size_t> edge_sizes;
    size_t begin_idx = 0;
    for (auto curr_size : idx_sizes) {
      size_t end_idx = begin_idx + curr_size;
      unsigned int begin_to_be_idx;
      f.seekg(header_size + begin_idx * sizeof(unsigned int));
      f.read(reinterpret_cast<char*>(&begin_to_be_idx), sizeof(begin_to_be_idx));
      unsigned int end_to_be_idx;
      f.seekg(header_size + end_idx * sizeof(unsigned int));
      f.read(reinterpret_cast<char*>(&end_to_be_idx), sizeof(end_to_be_idx));
      edge_sizes.push_back(end_to_be_idx - begin_to_be_idx);
      begin_idx = end_idx;
    }
    return edge_sizes;
  }

  struct read_partitioned_adj_part
    : hpx::parallel::detail::algorithm<read_partitioned_adj_part, int> {

    constexpr read_partitioned_adj_part() noexcept
      : hpx::parallel::detail::algorithm<read_partitioned_adj_part, int>(
          "read_partitioned_adj_part") {}

    template <typename ExPolicy>
    static int sequential(ExPolicy&&, partitioned_adjacency<0> G, size_t first_index,
                          size_t last_index, std::string file_name) {

      auto [indices, to_be_indexed] = deserialize_adj_part(file_name, first_index, last_index);
      if (last_index != G.size()) {
        indices.pop_back(); // Will be included in the next partition.
      }
      auto p_indices = G.get_indices().get_local_iterator(first_index).local();
      auto first_to_be_idx = indices.front();
      auto p_to_be_indexed =
        std::get<0>(G.get_to_be_indexed()).get_local_iterator(first_to_be_idx).local();

      std::copy(indices.begin(), indices.end(), p_indices);
      std::copy(to_be_indexed.begin(), to_be_indexed.end(), p_to_be_indexed);
      return 0;
    }

    template <typename ExPolicy>
    static int parallel(ExPolicy&&, partitioned_adjacency<0> G, size_t first_index,
                        size_t last_index, std::string file_name) {
      return 0;
    }
  };


  // struct read_partitioned_adj_part_action
  //   : hpx::actions::action<decltype(&read_partitioned_adj_part), &read_partitioned_adj_part,
  //                          read_partitioned_adj_part_action> {};


  auto partitioned_deserialize_adj(std::string mtx_file) {
    // Get sizes
    std::string file_name = detail::get_adj_filename(mtx_file);
    std::ifstream f(file_name, std::ios::binary | std::ios::in);
    char magic[30];
    f.read(magic, 30);
    if (strncmp(magic, "NWGRAPH ADJACENCY BINARY FILE", 30) != 0) {
      std::cerr << "Error: Invalid file format\n";
      exit(1);
    }

    size_t n_vertices;
    size_t n_edges;

    f.read(reinterpret_cast<char*>(&n_vertices), sizeof(n_vertices));
    f.read(reinterpret_cast<char*>(&n_edges), sizeof(n_edges));

    // construct partitioned adjacency graph
    size_t n_localities = hpx::get_num_localities(hpx::launch::sync);
    auto part_sizes = detail::partitioned_vertex_sizes(n_localities, n_vertices);
    auto edge_sizes = read_edge_sizes(file_name, part_sizes);
    partitioned_adjacency<0> G(n_vertices, n_edges, part_sizes, edge_sizes, "pg",
                               hpx::find_all_localities());

    partitioned_algorithm<read_partitioned_adj_part>(hpx::execution::seq, G, file_name);

    return G;
  }


} // namespace nw::graph
