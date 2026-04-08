#ifndef NW_GRAPH_PARTITIONED_SERIALIZE_HPP
#define NW_GRAPH_PARTITIONED_SERIALIZE_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>

#include "nwgraph/adjacency.hpp"
#include "nwgraph/distributed/algorithms/algorithm.hpp"
#include "nwgraph/edge_list.hpp"
#include "nwgraph/graph_base.hpp"
#include "nwgraph/io/mmio.hpp"
#include "nwgraph/distributed/adjacency.hpp"

#include <hpx/async_combinators/wait_all.hpp>
#include <hpx/executors/execution_policy.hpp>
#include <hpx/include/partitioned_vector_predef.hpp>
#include <hpx/parallel/segmented_algorithms/detail/dispatch.hpp>


namespace nw::graph {

  namespace detail {

    /* Returns a vector of partition sizes that sum to all_vertices.
       The sizes are as evenly distributed as possible.
    */
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

    /* Compresses a part of a edge list into a partial adjacency structure
     * This differs from a regular compression in that the vertex indices
     * are not assumed to start at 0.
     * Additionally, the vertex id cannot directly be used as an index into the
     * "indices" array, but must be offset by the minimum vertex id in the part.
     */
    auto compress_part(edge_list<directedness::undirected>& A, size_t idx_start, size_t n_vertices) {

      size_t n_edges = A.size();

      std::vector<default_index_t> indices(n_vertices, 0);
      std::vector<default_index_t> to_be_indexed(n_edges);


      // Need to sort edges by source vertex
      std::sort(A.begin(), A.end(), [](auto a, auto b) { return std::get<0>(a) < std::get<0>(b); });
      // Note that idx_start is not always the first vertex id in the part,
      // As some vertices may have no edges.

      size_t edge_idx = 0;
      for (auto&& [u, v] : A) {
        assert(u - idx_start < n_vertices);
        indices[u - idx_start]++;
        to_be_indexed[edge_idx++] = v - idx_start;
      }


      std::exclusive_scan(indices.begin(), indices.end(), indices.begin(), 0);

      return std::make_tuple(std::move(indices), std::move(to_be_indexed));
    };

    std::string get_adj_filename(std::string mtx_file) {
      return mtx_file.substr(0, mtx_file.find_last_of('.')) + ".adj";
    }


    template <typename id_t, typename vertex_id_t>
    class adj_writer {
      // The binary file format is as follows:
      // - A header containing :
      //      -A magic string "NWGRAPH ADJACENCY BINARY FILE"(30 bytes)
      //      -The number of vertices(size_t)
      //      - The number of edges(size_t)
      // - An array of indices(unsigned int) of size (n_vertices + 1)
      // - An array of to_be_indexed(unsigned int) of size (n_edges)
    public:
      adj_writer(std::string output_file, size_t n_vertices, size_t n_edges)
        : output_file_(output_file)
        , n_vertices_(n_vertices)
        , n_edges_(n_edges) {

        // Create file
        {
          std::ofstream f(output_file_, std::ios::binary);
        }

        // Open output file and initialize output stream positions
        auto open_mode = std::ios::binary | std::ios::in | std::ios::out | std::ios::ate;
        f_idx_.open(output_file_, open_mode);
        f_idx_.seekp(header_size_); // Move to the end of the header
        f_to_be_idx_.open(output_file_, open_mode);
        f_to_be_idx_.seekp(header_size_ + indices_size_); // Move to the end of the indices


        size_t total_size = header_size_ + indices_size_ + to_be_indexed_size_;
        std::filesystem::resize_file(output_file_,
                                     total_size); // Resize file to accommodate header and data
      }

    private:
      void write_header() {
        // We only write the header once we have written all data, to avoid
        // attempting to read incomplete files.
        std::fstream f(output_file_, std::ios::binary | std::ios::in | std::ios::out);
        f.write(magic_, sizeof(magic_));
        f.write(reinterpret_cast<const char*>(&n_vertices_), sizeof(n_vertices_));
        f.write(reinterpret_cast<const char*>(&n_edges_), sizeof(n_edges_));
        assert(f.tellp() == header_size_);
      }

    public:
      bool is_complete() {
        return (total_written_indices_ == n_vertices_ + 1) &&
          (total_written_to_be_indexed_ == n_edges_);
      }

      void write_next(std::vector<id_t>& indices, std::vector<vertex_id_t>& to_be_indexed) {
        auto f_idx_prev = f_idx_.tellp();
        // Serialize the indices and to_be_indexed vectors
        f_idx_.write(reinterpret_cast<const char*>(indices.data()), indices.size() * sizeof(id_t));

        assert((f_idx_.tellp() - f_idx_prev) / sizeof(id_t) == indices.size());
        assert(f_idx_.tellp() <= header_size_ + indices_size_);

        auto f_to_be_idx_prev = f_to_be_idx_.tellp();
        f_to_be_idx_.write(reinterpret_cast<const char*>(to_be_indexed.data()),
                           to_be_indexed.size() * sizeof(vertex_id_t));
        assert((f_to_be_idx_.tellp() - f_to_be_idx_prev) / sizeof(vertex_id_t) ==
               to_be_indexed.size());
        assert(f_to_be_idx_.tellp() <= header_size_ + indices_size_ + to_be_indexed_size_);

        total_written_indices_ += indices.size();
        total_written_to_be_indexed_ += to_be_indexed.size();

        if (is_complete()) {
          // Only write the header once all data has been written
          write_header();
        }
      }

    private:
      std::string output_file_;

      static constexpr char magic_[] = "NWGRAPH ADJACENCY BINARY FILE";

      size_t n_vertices_;
      size_t n_edges_;

      size_t header_size_ = sizeof(magic_) + sizeof(n_vertices_) + sizeof(n_edges_);
      size_t indices_size_ = sizeof(id_t) * (n_vertices_ + 1); // for adjacency list
      size_t to_be_indexed_size_ = sizeof(vertex_id_t) * (n_edges_); // for to_be_indexed

      // Keep counts for sanity checks
      size_t total_written_indices_ = 0;
      size_t total_written_to_be_indexed_ = 0;

      std::ofstream f_idx_;
      std::ofstream f_to_be_idx_;
    };


  } // namespace detail


  /* Serializes the adjacency of an undirected graph stored in a Matrix Market file into a binary
     (.adj) file. The parameter max_part_size allows processing graphs that do not fit in
     memory, loading and processing the graph in multiple partial passes.
  */
  std::string partitioned_serialize_adj(std::string mtx_file, size_t max_part_size = 2 << 26) {

    std::string file_name = detail::get_adj_filename(mtx_file);

    if (std::filesystem::exists(file_name)) {
      std::cout << "Adjacency file already exists, ";
      // See if magic number matches
      std::ifstream f(file_name);
      char magic[30];
      f.read(magic, sizeof(magic));
      if (strncmp(magic, "NWGRAPH ADJACENCY BINARY FILE", sizeof(magic)) == 0) {
        std::cout << "skipping serialization." << std::endl;
        return file_name;
      }
      std::cout << "but is invalid, overwriting." << std::endl;
    }

    // Get matrix size
    std::ifstream in_stream(mtx_file);
    auto [n_vertices, _, n_edges] = read_mm_metadata(in_stream);
    size_t n_parts = (n_vertices + max_part_size - 1) / max_part_size;
    auto part_sizes = detail::partitioned_vertex_sizes(n_parts, n_vertices);

    // Create adjacency writer
    using adj_writer_t = detail::adj_writer<default_index_t, default_vertex_id_type>;
    adj_writer_t writer(file_name, n_vertices, n_edges);

    size_t idx_start = 0;
    size_t offs_indices = 0; // Offset for indices
    for (size_t part_size : part_sizes) {

      in_stream = std::ifstream(mtx_file);

      // Read through the whole file, but only keeps edges whose source index is within the range
      size_t idx_end = idx_start + part_size;
      auto pred = [idx_start, idx_end](auto&& d0, auto&& d1, auto&& v)
      { return d0 >= idx_start && d0 < idx_end; };

      auto edgelist = read_mm<decltype(pred), directedness::undirected>(in_stream, pred);

      // Create partial adjacency
      auto [indices, to_be_indexed] = detail::compress_part(edgelist, idx_start, part_size);

      // Offset indices so that they are correct in the global context
      std::for_each(indices.begin(), indices.end(),
                    [offs_indices](auto& idx) { idx += offs_indices; });


      // Write parts to file
      writer.write_next(indices, to_be_indexed);

      // Update for next part
      offs_indices += to_be_indexed.size();
      idx_start = idx_end;
    }

    // Write final idx, which should point past the last element
    std::vector<default_index_t> final_idx = {(default_index_t)(n_edges)};
    std::vector<default_vertex_id_type> empty;
    writer.write_next(final_idx, empty);

    if (!writer.is_complete()) {
      std::cerr << "Error: Incomplete adjacency file written." << std::endl;
      exit(1);
    }

    size_t total_size = std::filesystem::file_size(file_name);
    std::cout << "Adjacency file serialized successfully." << std::endl;
    std::cout << "File size: " << total_size << " bytes." << std::endl;
    return file_name;
  }


  namespace detail {

    template <typename id_t, typename vertex_id_t>
    class adj_reader {
    public:
      adj_reader(std::string input_file)
        : input_file_(input_file) {
        // Read header
        std::ifstream f(input_file_, std::ios::binary | std::ios::in);
        char magic[30];
        f.read(magic, sizeof(magic));
        if (strncmp(magic, magic_, sizeof(magic_)) != 0) {
          std::cerr << "Error: Invalid file format\n";
          exit(1);
        }
        f.read(reinterpret_cast<char*>(&n_vertices_), sizeof(n_vertices_));
        f.read(reinterpret_cast<char*>(&n_edges_), sizeof(n_edges_));
        indices_size_ = sizeof(id_t) * (n_vertices_ + 1);
        to_be_indexed_size_ = sizeof(vertex_id_t) * (n_edges_);
      }

    public:
      std::pair<std::vector<id_t>, std::vector<vertex_id_t>> read_part(id_t begin_idx,
                                                                       id_t end_idx) {
        assert(begin_idx <= end_idx && end_idx <= n_vertices_ + 1);
        std::ifstream f(input_file_, std::ios::binary | std::ios::in);
        // Read indices
        std::vector<id_t> indices(end_idx - begin_idx + 1);
        f.seekg(header_size_ + begin_idx * sizeof(id_t));
        f.read(reinterpret_cast<char*>(indices.data()), indices.size() * sizeof(id_t));
        // Read to_be_indexed
        size_t begin_to_be_idx = indices.front();
        size_t end_to_be_idx = indices.back();
        std::vector<vertex_id_t> to_be_indexed(end_to_be_idx - begin_to_be_idx);
        f.seekg(header_size_ + indices_size_ + begin_to_be_idx * sizeof(vertex_id_t));
        f.read(reinterpret_cast<char*>(to_be_indexed.data()),
               to_be_indexed.size() * sizeof(vertex_id_t));

        return {std::move(indices), std::move(to_be_indexed)};
      }

      std::vector<size_t> get_edge_sizes(std::vector<size_t> idx_sizes) {
        std::ifstream f(input_file_, std::ios::binary | std::ios::in);
        f.seekg(header_size_ + indices_size_);
        std::vector<size_t> edge_sizes;
        size_t begin_idx = 0;
        for (auto curr_size : idx_sizes) {
          size_t end_idx = begin_idx + curr_size;
          id_t begin_to_be_idx;
          f.seekg(header_size_ + begin_idx * sizeof(id_t));
          f.read(reinterpret_cast<char*>(&begin_to_be_idx), sizeof(begin_to_be_idx));
          id_t end_to_be_idx;
          f.seekg(header_size_ + end_idx * sizeof(id_t));
          f.read(reinterpret_cast<char*>(&end_to_be_idx), sizeof(end_to_be_idx));
          edge_sizes.push_back(end_to_be_idx - begin_to_be_idx);
          begin_idx = end_idx;
        }
        return edge_sizes;
      }

      size_t n_vertices() const { return n_vertices_; }
      size_t n_edges() const { return n_edges_; }

    private:
      std::string input_file_;
      static constexpr char magic_[] = "NWGRAPH ADJACENCY BINARY FILE";
      size_t n_vertices_;
      size_t n_edges_;
      size_t header_size_ = sizeof(magic_) + sizeof(n_vertices_) + sizeof(n_edges_);
      size_t indices_size_;
      size_t to_be_indexed_size_;
    };


    struct read_partitioned_adj_part
      : hpx::parallel::detail::algorithm<read_partitioned_adj_part, int> {
      static constexpr bool partition_aware = true;

      constexpr read_partitioned_adj_part() noexcept
        : hpx::parallel::detail::algorithm<read_partitioned_adj_part, int>(
            "read_partitioned_adj_part") {}

      template <typename ExPolicy>
      static int sequential(ExPolicy&&, partitioned_adjacency<0> G,
                            const size_t first_index, const size_t last_index,
                            std::string file_name) {

        using reader_t = adj_reader<default_index_t, default_vertex_id_type>;
        reader_t reader(file_name);

        auto [indices, to_be_indexed] = reader.read_part(first_index, last_index);
        if (last_index < static_cast<std::size_t>(G.size())) {
          indices.pop_back();
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


  } // namespace detail


  auto partitioned_deserialize_adj(std::string mtx_file) {

    std::string file_name = detail::get_adj_filename(mtx_file);
    using reader_t = detail::adj_reader<default_index_t, default_vertex_id_type>;
    reader_t reader(file_name);

    size_t n_vertices = reader.n_vertices();
    size_t n_edges = reader.n_edges();

    // construct partitioned adjacency graph
    size_t n_localities = hpx::get_num_localities(hpx::launch::sync);
    auto part_sizes = detail::partitioned_vertex_sizes(n_localities, n_vertices);
    auto edge_sizes = reader.get_edge_sizes(part_sizes);

    partitioned_adjacency<0> G(n_vertices, n_edges, part_sizes, edge_sizes, "pg",
                               hpx::find_all_localities());

    partitioned_algorithm<detail::read_partitioned_adj_part>(hpx::execution::seq, G, file_name);

    return G;
  }


} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_SERIALIZE_HPP
