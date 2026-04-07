/**
 * @file partitioned_adjacency.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * @authors
 *   Andrew Lumsdaine
 *   Tony Liu
 *
 */

#ifndef NW_GRAPH_PARTITIONED_ADJACENCY_HPP
#define NW_GRAPH_PARTITIONED_ADJACENCY_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/adjacency.hpp"
#include "nwgraph/containers/partitioned_compressed.hpp"
#include "nwgraph/util/tag_invoke.hpp"

// #include "nwgraph/partitioned_build.hpp"

#include <algorithm>
#include <array>
#include <compare>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <vector>

#include <hpx/assert.hpp>
#include <hpx/include/partitioned_vector_predef.hpp>
#include <hpx/include/serialization.hpp>

namespace nw::graph {

  DECL_TAG_INVOKE(vertex_partition);
  DECL_TAG_INVOKE(partitions);
  DECL_TAG_INVOKE(local_view);
  DECL_TAG_INVOKE(parent);
  DECL_TAG_INVOKE(is_local_index);
  DECL_TAG_INVOKE(partition);
  DECL_TAG_INVOKE(remote_ref);

  class partition_descriptor;

  namespace detail {
    template <typename Graph>
    concept graph_with_index_partitions =
      requires(Graph const& G) {
        G.get_indices().partitions();
        G.get_indices().segment_begin();
        G.get_indices().segment_end();
      };

    hpx::id_type partition_id(partition_descriptor const& partition);
    partition_descriptor make_partition_descriptor(
      hpx::id_type partition_id, hpx::id_type locality, std::size_t first_index,
      std::size_t last_index);
  } // namespace detail

  class partition_descriptor {
    friend class hpx::serialization::access;

  public:
    partition_descriptor() = default;

    hpx::id_type locality() const { return locality_; }
    std::size_t first_index() const { return first_index_; }
    std::size_t last_index() const { return last_index_; }
    std::size_t size() const { return last_index_ - first_index_; }

    friend bool operator==(partition_descriptor lhs, partition_descriptor rhs) {
      return lhs.partition_id_ == rhs.partition_id_ && lhs.locality_ == rhs.locality_ &&
             lhs.first_index_ == rhs.first_index_ && lhs.last_index_ == rhs.last_index_;
    }

    friend auto operator<=>(partition_descriptor lhs, partition_descriptor rhs) {
      if (auto cmp = lhs.first_index_ <=> rhs.first_index_; cmp != 0) {
        return cmp;
      }
      if (auto cmp = lhs.last_index_ <=> rhs.last_index_; cmp != 0) {
        return cmp;
      }
      if (lhs.partition_id_ < rhs.partition_id_) {
        return std::strong_ordering::less;
      }
      if (rhs.partition_id_ < lhs.partition_id_) {
        return std::strong_ordering::greater;
      }
      if (lhs.locality_ < rhs.locality_) {
        return std::strong_ordering::less;
      }
      if (rhs.locality_ < lhs.locality_) {
        return std::strong_ordering::greater;
      }
      return std::strong_ordering::equal;
    }

  private:
    friend hpx::id_type detail::partition_id(partition_descriptor const&);
    friend partition_descriptor detail::make_partition_descriptor(
      hpx::id_type, hpx::id_type, std::size_t, std::size_t);

    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
      ar& partition_id_& locality_& first_index_& last_index_;
    }

    partition_descriptor(
      hpx::id_type partition_id, hpx::id_type locality, std::size_t first_index,
      std::size_t last_index)
      : partition_id_(HPX_MOVE(partition_id))
      , locality_(HPX_MOVE(locality))
      , first_index_(first_index)
      , last_index_(last_index) {}

    hpx::id_type partition_id_;
    hpx::id_type locality_;
    std::size_t first_index_ = 0;
    std::size_t last_index_ = 0;
  };

  namespace detail {
    inline hpx::id_type partition_id(partition_descriptor const& partition) {
      return partition.partition_id_;
    }

    template <typename PartitionedVector>
    auto find_partition_by_id(PartitionedVector const& pv, partition_descriptor partition) {
      return std::find_if(
        pv.segment_begin(), pv.segment_end(),
        [&](auto const& pv_partition) { return pv_partition.get_id() == partition_id(partition); });
    }

    template <typename PartitionedVector>
    auto find_matching_partition(PartitionedVector const& pv, partition_descriptor partition) {
      auto it = std::find_if(
        pv.segment_begin(), pv.segment_end(),
        [&](auto const& pv_partition)
        {
          auto first_index = static_cast<std::size_t>(pv_partition.first_);
          auto last_index = static_cast<std::size_t>(pv_partition.first_ + pv_partition.size_);
          last_index = std::min(last_index, static_cast<std::size_t>(pv.size()));
          return first_index == partition.first_index() && last_index == partition.last_index() &&
                 hpx::naming::get_locality_from_id(pv_partition.get_id()) == partition.locality();
        });
      HPX_ASSERT(it != pv.segment_end());
      return it;
    }

    template <typename PartitionedVector, typename SegmentIterator>
    partition_descriptor describe_partition(PartitionedVector const& pv, SegmentIterator partition) {
      auto partition_id = partition->get_id();
      auto locality = hpx::naming::get_locality_from_id(partition_id);
      auto first_index = static_cast<std::size_t>(partition->first_);
      auto last_index = static_cast<std::size_t>(partition->first_ + partition->size_);
      last_index = std::min(last_index, static_cast<std::size_t>(pv.size()));
      return make_partition_descriptor(
        HPX_MOVE(partition_id), HPX_MOVE(locality), first_index, last_index);
    }

    template <typename SourcePartitionedVector, typename TargetPartitionedVector>
    partition_descriptor aligned_partition(
      SourcePartitionedVector const& source, TargetPartitionedVector const& target,
      partition_descriptor partition) {
      auto source_partition = find_partition_by_id(source, partition);
      HPX_ASSERT(source_partition != source.segment_end());
      auto offset = std::distance(source.segment_begin(), source_partition);
      auto target_partition = target.segment_begin();
      std::advance(target_partition, offset);
      return describe_partition(target, target_partition);
    }

    inline partition_descriptor make_partition_descriptor(
      hpx::id_type partition_id, hpx::id_type locality, std::size_t first_index,
      std::size_t last_index) {
      return partition_descriptor(HPX_MOVE(partition_id), HPX_MOVE(locality), first_index,
                                  last_index);
    }
  }

  inline hpx::id_type partition_locality(partition_descriptor partition) {
    return partition.locality();
  }

#if 0
template <std::unsigned_integral index_type, std::unsigned_integral vertex_id_type, typename... Attributes>
class partitioned_index_compressed : public unipartite_graph_base, public partitioned_indexed_struct_of_arrays<index_type, vertex_id_type, Attributes...> {
  using base = partitioned_indexed_struct_of_arrays<index_type, vertex_id_type, Attributes...>;

public:
  using index_t = index_type;
  using vertex_id_type = vertex_id_type;

  partitioned_index_compressed(size_t N) : unipartite_graph_base(N), base(N) {}

  void close_for_push_back() { base::close_for_push_back(); };

  auto num_edges() { return base::to_be_indexed_.size(); }
};

template <typename... Attributes>
using compressed = partitioned_index_compressed<default_index_t, default_vertex_id_type, Attributes...>;
#endif

  /**
   * @brief Index partitioned_adjacency structure. This data structures stores unipartite graph in
   * Compressed Sparse Row format. The underlying data structure is a structure of arrays for
   * storage. Index_adjacency is a partitioned_adjacency list representation of a unipartite
   * graph. A unipartite graph has one vertex set.
   *
   * @tparam idx The index type to indicate the type of the graph, can be either 0 or 1.
   * @tparam index_type The data type used to represent a vertex index, required to be os unsigned
   * integral type.
   * @tparam vertex_id The data type used to represent a vertex ID, required to be os unsigned
   * integral type.
   * @tparam Attributes A variadic list of edge property types.
   */
  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
            typename... Attributes>
  class partitioned_index_adjacency
    : public unipartite_graph_base,
      public partitioned_indexed_struct_of_arrays<index_type, vertex_id, Attributes...> {
    using base = partitioned_indexed_struct_of_arrays<index_type, vertex_id, Attributes...>;

    friend class hpx::serialization::access;

    // clang-format off
    void serialize(hpx::serialization::input_archive& ar, unsigned) {
      ar >> hpx::serialization::base_object<unipartite_graph_base>(*this)
         >> hpx::serialization::base_object<base>(*this);
    }

    void serialize(hpx::serialization::output_archive& ar, unsigned) const {
      ar << hpx::serialization::base_object<unipartite_graph_base>(*this)
         << hpx::serialization::base_object<base>(*this);
    }
    // clang-format on

    // shallow copy constructor, shallow-copies partitioned vectors
    partitioned_index_adjacency(partitioned_index_adjacency const& rhs, bool make_unmanaged)
      : unipartite_graph_base(static_cast<unipartite_graph_base const&>(rhs))
      , base(static_cast<base const&>(rhs), make_unmanaged) {}

    template <std::size_t I>
    static std::string generate_name(std::string const& name) {
      return name + std::to_string(I);
    }

  public:
    using index_t = index_type;
    using vertex_id_type = vertex_id;
    using num_vertices_type = std::array<size_t, 1>;
    using num_edges_type = index_t;

    // The first index_t isn't considered an attribute.
    using attributes_t = std::tuple<Attributes...>;
    static constexpr std::size_t getNAttr() { return sizeof...(Attributes); }

    // needed for HPX serialization
    partitioned_index_adjacency() = default;

    /**
     * @brief Constructor of partitioned_index_adjacency. Require the type of the graph to be
     * unipartite. Create an empty partitioned_index_adjacency.
     */
    template <typename Vector>
    partitioned_index_adjacency(
      size_t N, size_t M, Vector&& index_sizes, Vector&& to_be_index_sizes, char const* name = "pg",
      std::vector<hpx::id_type> const& localities = hpx::find_all_localities())
      requires(std::is_same_v<unipartite_graph_base, unipartite_graph_base>)
      : unipartite_graph_base(N)
      , base(N, N + 1, M, std::forward<Vector>(index_sizes),
             std::forward<Vector>(to_be_index_sizes), generate_name<idx>(name), localities) {}

    /**
     * @brief Constructor of partitioned_index_adjacency. Require the type of the graph to be
     * unipartite. Copies loc_adj to the partitioned_index_adjacency.
     */
    template <typename Vector>
    partitioned_index_adjacency(
      size_t N, size_t M, Vector&& index_sizes, Vector&& to_be_index_sizes,
      adjacency<idx>& local_adj, char const* name = "pg",
      std::vector<hpx::id_type> const& localities = hpx::find_all_localities())
      requires(std::is_same_v<unipartite_graph_base, unipartite_graph_base>)
      : unipartite_graph_base(N)
      , base(N, N + 1, M, std::forward<Vector>(index_sizes),
             std::forward<Vector>(to_be_index_sizes), local_adj, generate_name<idx>(name),
             localities) {}

    template <typename Vector>
    partitioned_index_adjacency(
      size_t N, size_t N1, size_t M, Vector&& index_sizes, Vector&& to_be_index_sizes,
      char const* name = "pg",
      std::vector<hpx::id_type> const& localities = hpx::find_all_localities())
      requires(std::is_same_v<unipartite_graph_base, unipartite_graph_base>)
      : unipartite_graph_base(N)
      , base(N, N1, M, std::forward<Vector>(index_sizes), std::forward<Vector>(to_be_index_sizes),
             generate_name<idx>(name), localities) {}

    // Create reference to partitioned_index_adjacency
    partitioned_index_adjacency ref() const { return partitioned_index_adjacency(*this, true); }

    /**
     * @brief Constructor of partitioned_index_adjacency. Require the type of the graph to be
     * unipartite. Create an empty partitioned_index_adjacency.
     */
    //  partitioned_index_adjacency(std::array<size_t, 1> N, size_t M = 0)
    //  requires(std::is_same<unipartite_graph_base, unipartite_graph_base>::value) :
    //  unipartite_graph_base(N), base(N[0], M) {}

    template <typename Vector, class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
    partitioned_index_adjacency(
      index_edge_list<vertex_id_type, unipartite_graph_base, directedness::directed, Attributes...>&
        A,
      Vector&& index_sizes, Vector&& to_be_index_sizes, bool sort_adjacency = false,
      ExecutionPolicy&& policy = {}, char const* name = "pg",
      std::vector<hpx::id_type> const& localities = hpx::find_all_localities())
      requires(std::is_same_v<unipartite_graph_base, unipartite_graph_base>)
      : unipartite_graph_base(A.num_vertices()[0])
      , base(A.num_vertices()[0] + 1, A.num_vertices()[0] + 1, A.num_edges(),
             std::forward<Vector>(index_sizes), std::forward<Vector>(to_be_index_sizes),
             generate_name<idx>(name), localities) {
      // partitioned_fill<idx>(A, *this, sort_adjacency, policy);
    }

    // template <class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
    // partitioned_index_adjacency(index_edge_list<vertex_id_type, unipartite_graph_base,
    //                                             directedness::undirected, Attributes...>& A,
    //                             bool sort_adjacency = false, ExecutionPolicy&& policy = {})
    //   requires(std::is_same_v<unipartite_graph_base, unipartite_graph_base>)
    //   : unipartite_graph_base(A.num_vertices()[0])
    //   , base(A.num_vertices()[0] + 1) {
    //   fill<idx>(A, *this, sort_adjacency, policy);
    // }

    // template <class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
    // partitioned_index_adjacency(size_t N,
    //                             index_edge_list<vertex_id_type, unipartite_graph_base,
    //                                             directedness::directed, Attributes...>& A,
    //                             bool sort_adjacency = false, ExecutionPolicy&& policy = {})
    //   requires(std::is_same_v<unipartite_graph_base, unipartite_graph_base>)
    //   : unipartite_graph_base(N)
    //   , base(N) {
    //   fill<idx>(A, *this, sort_adjacency, policy);
    // }
    // template <class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
    // partitioned_index_adjacency(size_t N,
    //                             index_edge_list<vertex_id_type, unipartite_graph_base,
    //                                             directedness::undirected, Attributes...>& A,
    //                             bool sort_adjacency = false, ExecutionPolicy&& policy = {})
    //   requires(std::is_same_v<unipartite_graph_base, unipartite_graph_base>)
    //   : unipartite_graph_base(N)
    //   , base(N) {
    //   fill<idx>(A, *this, sort_adjacency, policy);
    // }
    //// customized move constructor
    // partitioned_index_adjacency(std::vector<vertex_id>&& indices,
    //                             std::vector<vertex_id>&& first_to_be,
    //                             std::vector<Attributes>&&... rest_to_be)
    //   requires(std::is_same_v<unipartite_graph_base, unipartite_graph_base>)
    //   : unipartite_graph_base(indices.size() - 1)
    //   , base(std::move(indices), std::move(first_to_be), std::move(rest_to_be)...) {}
    // partitioned_index_adjacency(
    //   std::vector<vertex_id>&& indices,
    //   std::tuple<std::vector<vertex_id>, std::vector<Attributes>...>&& to_be_indexed)
    //   requires(std::is_same<unipartite_graph_base, unipartite_graph_base>::value)
    //   : unipartite_graph_base(indices.size() - 1)
    //   , base(std::move(indices), std::move(to_be_indexed)) {}
    //// customized copy constructor
    // partitioned_index_adjacency(const std::vector<vertex_id>& indices,
    //                             const std::vector<vertex_id>& first_to_be,
    //                             const std::vector<Attributes>&... rest_to_be)
    //   requires(std::is_same<unipartite_graph_base, unipartite_graph_base>::value)
    //   : unipartite_graph_base(indices.size() - 1)
    //   , base(indices, first_to_be, rest_to_be...) {}
    // partitioned_index_adjacency(
    //   const std::vector<vertex_id>& indices,
    //   const std::tuple<std::vector<vertex_id>, std::vector<Attributes>...>& to_be_indexed)
    //   requires(std::is_same<unipartite_graph_base, unipartite_graph_base>::value)
    //   : unipartite_graph_base(indices.size() - 1)
    //   , base(indices, to_be_indexed) {}

    num_vertices_type num_vertices() const { return {base::size()}; };
    num_edges_type num_edges() const { return base::to_be_indexed_.size(); };
    /**
     * @brief Serialize the partitioned_index_adjacency into binary file.
     *
     * @param outfile_name The output file name.
     */
    void serialize(const std::string& outfile_name) const {
      std::ofstream out_file(outfile_name, std::ofstream::binary);
      unipartite_graph_base::serialize(out_file);
      base::serialize(out_file);
    }

    /**
     * @brief Deserialize the binary into partitioned_index_adjacency.
     *
     * @param infile_name The input file name.
     */
    void deserialize(const std::string& infile_name) {
      std::ifstream infile(infile_name, std::ifstream::binary);
      unipartite_graph_base::deserialize(infile);
      base::deserialize(infile);
    }
  };

  template <int idx, typename... Attributes>
  using partitioned_adjacency =
    partitioned_index_adjacency<idx, default_index_t, default_vertex_id_type, Attributes...>;

  namespace detail {
    template <typename Graph, typename SegmentIterator>
    partition_descriptor partition_from_segment(Graph const& G, SegmentIterator partition) {
      return describe_partition(G, partition);
    }
  } // namespace detail

  template <typename Graph, std::unsigned_integral VertexId>
    requires detail::graph_with_index_partitions<std::remove_reference_t<Graph>> &&
             std::convertible_to<VertexId, vertex_id_t<std::remove_reference_t<Graph>>>
  partition_descriptor tag_invoke(vertex_partition_tag, Graph const& G, VertexId v) {
    auto partition = G.get_indices().get_segment_iterator(static_cast<std::size_t>(v));
    return detail::partition_from_segment(G, partition);
  }

  template <typename Graph>
    requires detail::graph_with_index_partitions<std::remove_reference_t<Graph>>
  std::vector<partition_descriptor> tag_invoke(partitions_tag, Graph const& G) {
    std::vector<partition_descriptor> result;
    auto const& graph_partitions = G.get_indices().partitions();
    result.reserve(graph_partitions.size());

    auto partition = G.get_indices().segment_begin();
    auto partition_end = G.get_indices().segment_end();
    for (; partition != partition_end; ++partition) {
      result.push_back(detail::partition_from_segment(G, partition));
    }
    return result;
  }

  template <std::unsigned_integral VertexId>
  bool is_local(partition_descriptor partition, VertexId v) {
    auto index = static_cast<std::size_t>(v);
    return partition.first_index() <= index && index < partition.last_index();
  }

  template <typename G>
  concept partitioned_graph =
    requires(std::remove_reference_t<G> const& g,
             vertex_id_t<std::remove_reference_t<G>> v,
             partition_descriptor partition) {
      typename vertex_id_t<std::remove_reference_t<G>>;
      { g.size() } -> std::convertible_to<std::size_t>;
      { vertex_partition(g, v) } -> std::same_as<partition_descriptor>;
      { partitions(g) } -> std::ranges::random_access_range;
      requires std::same_as<std::ranges::range_value_t<decltype(partitions(g))>,
                            partition_descriptor>;
      { partition_locality(partition) } -> std::same_as<hpx::id_type>;
      { is_local(partition, v) } -> std::convertible_to<bool>;
    };

  template <typename G>
  concept copartitioned_graph =
    partitioned_graph<G> &&
    requires(std::remove_reference_t<G> const& g) {
      { g.get_indices().size() } -> std::convertible_to<std::size_t>;
      { g.get_indices().get_partition_sizes() };
      { g.get_indices().get_partition_localities() };
    };

  template <int idx, edge_list_graph edge_list_t>
  auto make_partitioned_adjacency(edge_list_t& el) {
    return partitioned_adjacency<idx>(el);
  }


  template <int idx, edge_list_c edge_list_t, std::unsigned_integral u_integral,
            class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
  auto make_partitioned_adjacency(edge_list_t& el, u_integral n,
                                  directedness edge_directedness = directedness::directed,
                                  ExecutionPolicy&& policy = {}) {
    partitioned_adjacency<idx> adj(n);
    fill<idx>(el, adj, edge_directedness, policy);
    return adj;
  }

  /**
   * @brief Index biadjacency structure. This data structures stores bipartite graph in Compressed
   * Sparse Row format. The underlying data structure is a structure of arrays for storage.
   * Index_biadjacency is an biadjacency list represenation of a bipartite graph.
   * A bipartite graph has two vertex partitions.
   *
   *
   * @tparam idx The index type to indicate the type of the graph, can be either 0 or 1.
   * @tparam index_type The data type used to represent a vertex index, required to be os unsigned
   * integral type.
   * @tparam vertex_id The data type used to represent a vertex ID, required to be os unsigned
   * integral type.
   * @tparam Attributes A variadic list of edge property types.
   */
  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
            typename... Attributes>
  class partitioned_index_biadjacency
    : public bipartite_graph_base,
      public hpx::partitioned_vector<
        index_type, indexed_struct_of_arrays<index_type, vertex_id, Attributes...>> {
    using base =
      hpx::partitioned_vector<index_type,
                              indexed_struct_of_arrays<index_type, vertex_id, Attributes...>>;

  public:
    using value_type = index_type;
    using index_t = index_type;
    using vertex_id_type = vertex_id;
    using num_vertices_type = std::array<vertex_id_type, 2>;
    using num_edges_type = index_t;

    // The first index_t isn't considered an attribute.
    using attributes_t = std::tuple<Attributes...>;
    static constexpr std::size_t getNAttr() { return sizeof...(Attributes); }

    /**
     * @brief Constructor of index_biadjacency. Require the type of the graph to be bipartite.
     * Create an empty index_biadjacency.
     */
    partitioned_index_biadjacency(size_t N0 = 0, size_t N1 = 0, size_t M = 0)
      requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
      : bipartite_graph_base(N0, N1)
      , base(N0, M) {}
    /**
     * @brief Constructor of index_biadjacency. Require the type of the graph to be bipartite.
     * Create an empty index_biadjacency.
     */
    partitioned_index_biadjacency(std::array<size_t, 2> N, size_t M = 0)
      requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
      : bipartite_graph_base(N[idx], N[(idx + 1) % 2])
      , base(N[idx], M) {}

    template <class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
    partitioned_index_biadjacency(index_edge_list<vertex_id_type, bipartite_graph_base,
                                                  directedness::directed, Attributes...>& A,
                                  bool sort_biadjacency = false, ExecutionPolicy&& policy = {})
      requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
      : bipartite_graph_base(A.num_vertices()[idx], A.num_vertices()[(idx + 1) % 2])
      , base(A.num_vertices()[idx] + 1) {
      fill_biadjacency<idx>(A, *this, sort_biadjacency, policy);
    }

    template <class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
    partitioned_index_biadjacency(index_edge_list<vertex_id_type, bipartite_graph_base,
                                                  directedness::undirected, Attributes...>& A,
                                  bool sort_biadjacency = false, ExecutionPolicy&& policy = {})
      requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
      : bipartite_graph_base(A.num_vertices()[idx], A.num_vertices()[(idx + 1) % 2])
      , base(A.num_vertices()[idx] + 1) {
      fill_biadjacency<idx>(A, *this, sort_biadjacency, policy);
    }
    // customized move constructor
    partitioned_index_biadjacency(size_t N1, std::vector<vertex_id>&& indices,
                                  std::vector<vertex_id>&& first_to_be,
                                  std::vector<Attributes>&&... rest_to_be)
      requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
      : bipartite_graph_base(indices.size() - 1, N1)
      , base(std::move(indices), std::move(first_to_be), std::move(rest_to_be)...) {}
    partitioned_index_biadjacency(
      size_t N1, std::vector<vertex_id>&& indices,
      std::tuple<std::vector<vertex_id>, std::vector<Attributes>...>&& to_be_indexed)
      : bipartite_graph_base(indices.size() - 1, N1)
      , base(std::move(indices), std::move(to_be_indexed)) {}
    // customized copy constructor
    partitioned_index_biadjacency(size_t N1, const std::vector<vertex_id>& indices,
                                  const std::vector<vertex_id>& first_to_be,
                                  const std::vector<Attributes>&... rest_to_be)
      requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
      : bipartite_graph_base(indices.size() - 1, N1)
      , base(indices, first_to_be, rest_to_be...) {}
    partitioned_index_biadjacency(
      size_t N1, const std::vector<vertex_id>& indices,
      const std::tuple<std::vector<vertex_id>, std::vector<Attributes>...>& to_be_indexed)
      requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
      : bipartite_graph_base(indices.size() - 1, N1)
      , base(indices, to_be_indexed) {}

    auto num_vertices() const { return vertex_cardinality; }
    num_edges_type num_edges() const { return base::to_be_indexed_.size(); };
    /**
     * @brief Serialize the partitioned_index_adjacency into binary file.
     *
     * @param outfile_name The output file name.
     */
    void serialize(const std::string& outfile_name) {
      std::ofstream out_file(outfile_name, std::ofstream::binary);
      bipartite_graph_base::serialize(out_file);
      base::serialize(out_file);
    }

    /**
     * @brief Deserialize the binary into partitioned_index_adjacency.
     *
     * @param infile_name The input file name.
     */
    void deserialize(const std::string& infile_name) {
      std::ifstream infile(infile_name, std::ifstream::binary);
      bipartite_graph_base::deserialize(infile);
      base::deserialize(infile);
    }
  };

  template <int idx, typename... Attributes>
  using partitioned_biadjacency =
    partitioned_index_biadjacency<idx, default_index_t, default_vertex_id_type, Attributes...>;

  template <int idx, edge_list_graph edge_list_t>
  auto make_partitioned_biadjacency(edge_list_t& el) {
    return partitioned_biadjacency<idx>(el);
  }


  template <int idx, edge_list_c edge_list_t, std::unsigned_integral u_integral,
            class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
  auto make_partitioned_biadjacency(edge_list_t& el, u_integral n0, u_integral n1,
                                    directedness edge_directedness = directedness::directed,
                                    ExecutionPolicy&& policy = {}) {
    biadjacency<idx> adj(n0, n1);
    fill_biadjacency<idx>(el, adj, policy);
    return adj;
  }


  // template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
  // typename... Attributes> auto num_vertices(const partitioned_index_adjacency<idx, index_type,
  // vertex_id_type, Attributes...>& g) {
  //   return g.num_vertices();
  // }
  // partitioned_index_adjacency num_vertices CPO
  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
            typename... Attributes>
  auto
  tag_invoke(const num_vertices_tag,
             const partitioned_index_adjacency<idx, index_type, vertex_id_type, Attributes...>& g) {
    return g.num_vertices()[0];
  }
  // partitioned_index_adjacency degree CPO
  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
            std::unsigned_integral lookup_type, typename... Attributes>
  auto
  tag_invoke(const degree_tag,
             const partitioned_index_adjacency<idx, index_type, vertex_id_type, Attributes...>& g,
             lookup_type i) {
    return g[i].size();
  }
  // partitioned_index_adjacency degree CPO
  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
            typename... Attributes>
  auto
  tag_invoke(const degree_tag,
             const partitioned_index_adjacency<idx, index_type, vertex_id_type, Attributes...>& g,
             const typename partitioned_index_adjacency<idx, index_type, vertex_id_type,
                                                        Attributes...>::sub_view& v) {
    return v.size();
  }
  // partitioned_index_biadjacency num_vertices CPO
  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
            typename... Attributes>
  auto
  tag_invoke(const num_vertices_tag,
             const partitioned_index_biadjacency<idx, index_type, vertex_id_type, Attributes...>& g,
             int jdx = 0) {
    return g.num_vertices()[jdx];
  }
  // partitioned_index_biadjacency degree CPO
  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
            std::unsigned_integral lookup_type, typename... Attributes>
  auto
  tag_invoke(const degree_tag,
             const partitioned_index_biadjacency<idx, index_type, vertex_id_type, Attributes...>& g,
             lookup_type i) {
    return g[i].size();
  }
  // partitioned_index_biadjacency degree CPO
  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
            typename... Attributes>
  auto
  tag_invoke(const degree_tag,
             const partitioned_index_biadjacency<idx, index_type, vertex_id_type, Attributes...>& g,
             const typename partitioned_index_biadjacency<idx, index_type, vertex_id_type,
                                                          Attributes...>::sub_view& v) {
    return v.size();
  }

} // namespace nw::graph

// serialization of partitioned_index_adjacency requires special handling
template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
          typename... Attributes>
struct hpx::traits::needs_reference_semantics<
  nw::graph::partitioned_index_adjacency<idx, index_type, vertex_id, Attributes...>>
  : std::true_type {};

#endif // NW_GRAPH_PARTITIONED_ADJACENCY_HPP
