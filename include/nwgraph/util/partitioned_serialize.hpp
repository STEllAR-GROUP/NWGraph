#include <filesystem>

#include "nwgraph/adjacency.hpp"
#include "nwgraph/edge_list.hpp"
#include "nwgraph/graph_base.hpp"
#include "nwgraph/io/mmio.hpp"
#include "nwgraph/partitioned_adjacency.hpp"

//using unsigned_int = unsigned int;
//HPX_REGISTER_PARTITIONED_VECTOR(unsigned_int)

namespace nw::graph {

  std::tuple<size_t, size_t> get_part_range(size_t nNonzeros, size_t n_partitions, size_t i) {
    size_t segment_size = (nNonzeros + n_partitions - 1) / n_partitions;
    size_t begin = i * segment_size;
    size_t end = std::min((i + 1) * segment_size, nNonzeros);
    return {begin, end};
  }

  std::string get_part_filename(std::string in_file, size_t n_partitions, size_t i) {
    return in_file + "." + std::to_string(i) + "." + std::to_string(n_partitions) + ".bmtk";
  }


  std::string get_index_filename(std::string in_file, size_t n_partitions) {
    return in_file + ".index." + std::to_string(n_partitions) + ".bmtk";
  }


  bool segment_files_exist(std::string in_file, size_t n_partitions) {
    for (size_t i = 0; i < n_partitions; ++i) {
      std::string out_file = get_part_filename(in_file, n_partitions, i);
      if (!std::filesystem::exists(out_file)) {
        return false;
      }
    }
    std::string index_file = get_index_filename(in_file, n_partitions);
    if (!std::filesystem::exists(index_file)) {
      return false;
    }
    return true;
  }


  // Serialize vector
  template <typename T>
  std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
    out << v.size() << "\n";
    for (const auto& e : v) {
      out << e << "\n";
    }
    return out;
  }


  // Deserialize vector
  template <typename T>
  std::istream& operator>>(std::istream& in, std::vector<T>& v) {
    size_t size;
    in >> size;
    v.resize(size);
    for (auto& e : v) {
      in >> e;
    }
    return in;
  }


  auto compress_part(edge_list<directedness::undirected>& A, size_t n_indices,
                     size_t n_to_be_indexed, size_t prev_index_end, bool is_last) {
   
    std::vector<default_index_t> indices(n_indices + 1, 0);
    std::vector<default_index_t> to_be_indexed(n_to_be_indexed);

    size_t idx_offset = std::get<0>(A[0]);
    for (size_t i = 0; i < A.size(); i++) {
      indices[std::get<0>(A[i]) - idx_offset]++;
      to_be_indexed[i] = std::get<1>(A[i]);
    }

    std::exclusive_scan(indices.begin(), indices.end(), indices.begin(), prev_index_end);

    // Last element only needed for the last segment
    if (!is_last) {
      indices.pop_back();
    }

    return std::make_tuple(std::move(indices), std::move(to_be_indexed));
  };


  void serialize_adjacency_graph_segments(std::string mtx_file, size_t n_partitions) {

    // Get matrix size
    std::ifstream in_stream(mtx_file);
    auto [n_vertices, _, n_edges] = read_mm_metadata(in_stream);

    // Keep track of the number of vertices and edges in each segment
    std::vector<size_t> vertex_counts;
    std::vector<size_t> edge_counts;

    size_t prev_index_end = 0;

    for (size_t i = 0; i < n_partitions; ++i) {
      auto [index_begin, index_end] = get_part_range(n_vertices, n_partitions, i);

      in_stream = std::ifstream(mtx_file);

      // Reads the whole file, but only keeps edges whose source index is within the range
      auto pred = [index_begin, index_end](auto&& d0, auto&& d1, auto&& v)
      { return d0 >= index_begin && d0 < index_end; };

      auto edgelist = read_mm<decltype(pred), directedness::undirected>(in_stream, pred);

      size_t n_indices = index_end - index_begin;
      size_t n_to_be_indexed = edgelist.size();

      vertex_counts.push_back(n_indices);
      edge_counts.push_back(n_to_be_indexed);

      // compress the segment
      sort_by<0>(edgelist);
      bool is_last = i == n_partitions - 1;
      auto [indices, to_be_indexed] =
        compress_part(edgelist, n_indices, n_to_be_indexed, prev_index_end, is_last);

      prev_index_end += n_to_be_indexed;

      // Serialize the compressed segment
      std::string out_file = get_part_filename(mtx_file, n_partitions, i);

      std::ofstream out_stream(out_file, std::ofstream::binary);
      out_stream << mtx_file << "\n";
      out_stream << index_begin << "\n";
      out_stream << index_end << "\n";
      out_stream << n_indices << "\n";
      out_stream << n_to_be_indexed << "\n";
      out_stream << i << "\n";
      out_stream << n_partitions << "\n";
      out_stream << indices;
      out_stream << to_be_indexed;
    }

    // Write the index file
    // Index contains the number of vertices and edges in each segment, needed
    // to initialize the partitioned adjacency graph
    std::string index_filename = get_index_filename(mtx_file, n_partitions);
    std::ofstream index_stream(index_filename);
    index_stream << mtx_file << "\n";
    index_stream << n_vertices << "\n";
    index_stream << n_edges << "\n";
    index_stream << n_partitions << "\n";
    index_stream << vertex_counts << edge_counts;
  }


  auto deserialize_adjacency_graph_segment(std::string mtx_file, size_t n_partitions, size_t idx) {
    std::string out_file = get_part_filename(mtx_file, n_partitions, idx);

    std::ifstream in_stream(out_file, std::ifstream::binary);
    std::string origin_file;
    size_t segment_begin, segment_end, n_indices, n_to_be_indexed, i_segment_idx, i_n_partitions;
    in_stream >> origin_file >> segment_begin >> segment_end >> n_indices >> n_to_be_indexed >>
      i_segment_idx >> i_n_partitions;

    if (origin_file != mtx_file) {
      std::cerr << "Error: origin file mismatch\n";
      exit(1);
    }

    if (idx != i_segment_idx) {
      std::cerr << "Error: segment index mismatch\n";
      exit(1);
    }

    if (n_partitions != i_n_partitions) {
      std::cerr << "Error: segment count mismatch\n";
      exit(1);
    }
    
    std::vector<default_index_t> indices;
    std::vector<default_index_t> to_be_indexed;

    in_stream >> indices;
    in_stream >> to_be_indexed;

    return std::make_tuple(std::move(indices), std::move(to_be_indexed));
  }


  auto construct_partitioned_adjacency_from_index(std::string index_file, size_t n_partitions) {

    // Read sizes from index file
    std::vector<size_t> vert_sizes(n_partitions);
    std::vector<size_t> edge_sizes(n_partitions);

    if (!std::filesystem::exists(index_file)) {
      std::cerr << "Error: index file (" + index_file + ") does not exist\n ";
      exit(1);
    }

    std::ifstream index_stream(index_file);
    std::string origin_file;
    size_t n_vertices, n_edges, i_n_partitions;
    index_stream >> origin_file >> n_vertices >> n_edges >> i_n_partitions >> vert_sizes >>
      edge_sizes;

    // Construct the partitioned adjacency graph
    partitioned_adjacency<0> G(n_vertices, n_edges, vert_sizes, edge_sizes, "pg",
                               hpx::find_all_localities());

    return G;
  }

  partitioned_adjacency<0> load_partitioned_adjacency(std::string mtx_file, size_t n_partitions) {

    std::string index_file = get_index_filename(mtx_file, n_partitions);
    auto G = construct_partitioned_adjacency_from_index(index_file, n_partitions);

    auto& p_indices = G.get_indices();
    auto& p_to_be_indexed = G.get_to_be_indexed();

    // Fill partitioned adjacency graph
    auto it_indices = p_indices.begin();
    auto it_to_be_indexed = p_to_be_indexed.begin();
    for (size_t i = 0; i < n_partitions; i++) {

      std::string filename = get_part_filename(mtx_file, n_partitions, i);

      // Deserialize the segment
      auto [indices, to_be_indexed] =
        deserialize_adjacency_graph_segment(mtx_file, n_partitions, i);

      // TODO: Assert that partition sizes match

      std::vector<std::size_t> const empty;
      p_indices.set_values(hpx::launch::sync, i, empty, HPX_MOVE(indices));

      std::get<0>(p_to_be_indexed).set_values(hpx::launch::sync, i, empty, HPX_MOVE(to_be_indexed));


    }

    return G;
  }
} // namespace nw::graph
