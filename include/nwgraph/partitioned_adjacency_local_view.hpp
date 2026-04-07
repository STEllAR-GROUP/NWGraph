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

#ifndef NW_GRAPH_PARTITIONED_ADJACENCY_LOCAL_VIEW_HPP
#define NW_GRAPH_PARTITIONED_ADJACENCY_LOCAL_VIEW_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/partitioned_adjacency.hpp"
#include "nwgraph/containers/partitioned_compressed_local_view.hpp"
#include "nwgraph/containers/partitioned_soa_local_view.hpp"
#include "nwgraph/util/tag_invoke.hpp"

// #include "nwgraph/partitioned_build.hpp"

#include <array>
#include <cstddef>
#include <concepts>
#include <type_traits>

#include <hpx/assert.hpp>
#include <hpx/include/partitioned_vector_predef.hpp>

namespace nw::graph {

  namespace detail {
    template <typename Graph>
    struct local_graph_view;
  } // namespace detail

  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
            typename... Attributes>
  class partitioned_index_adjacency_local_view
    : public unipartite_graph_base, // Haven't checked if all properties are satisfied, but this
                                    // allows us to use some useful CPOs (num_vertices, target, etc.)
      public partitioned_indexed_struct_of_arrays_local_view<index_type, vertex_id, Attributes...> {
    using base = partitioned_indexed_struct_of_arrays_local_view<index_type, vertex_id, Attributes...>;

    partitioned_index_adjacency<idx, index_type, vertex_id, Attributes...>* parent_ = nullptr;
  public:
    using index_t = index_type;
    using vertex_id_type = vertex_id;
    using num_vertices_type = size_t;
    using num_edges_type = index_t;

    using graph_type = partitioned_index_adjacency<idx, index_type, vertex_id, Attributes...>;

    // The first index_t isn't considered an attribute.
    using attributes_t = std::tuple<Attributes...>;
    static constexpr std::size_t getNAttr() { return sizeof...(Attributes); }

    partitioned_index_adjacency_local_view(
      partitioned_index_adjacency<idx, index_type, vertex_id, Attributes...>& parent,
      partition_descriptor partition)
      : base(parent, partition)
      , parent_(&parent) {}


    num_vertices_type num_vertices() const { return {base::size()}; };
    num_edges_type num_edges() const { return base::to_be_indexed_.size(); };

    graph_type& parent() { 
        assert(parent_ != nullptr);
        return *parent_; 
    }

    graph_type const& parent() const {
        assert(parent_ != nullptr);
        return *parent_;
    }
  };

  template <int idx, typename... Attributes>
  using partitioned_adjacency_local_view =
    partitioned_index_adjacency_local_view<idx, default_index_t, default_vertex_id_type,
                                           Attributes...>;

  namespace detail {
    template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
              typename... Attributes>
    struct local_graph_view<partitioned_index_adjacency<idx, index_type, vertex_id, Attributes...>> {
      using type =
        partitioned_index_adjacency_local_view<idx, index_type, vertex_id, Attributes...>;
    };

  } // namespace detail

  template <typename Graph>
    requires requires { typename detail::local_graph_view<std::remove_reference_t<Graph>>::type; }
  auto tag_invoke(local_view_tag, Graph& G, partition_descriptor partition)
    -> typename detail::local_graph_view<std::remove_reference_t<Graph>>::type {
    using local_view_t = typename detail::local_graph_view<std::remove_reference_t<Graph>>::type;
    return local_view_t(G, partition);
  }

  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
            typename... Attributes>
  auto& tag_invoke(parent_tag,
    partitioned_index_adjacency_local_view<idx, index_type, vertex_id, Attributes...>& G) {
    return G.parent();
  }

  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
            typename... Attributes>
  auto const& tag_invoke(parent_tag,
    partitioned_index_adjacency_local_view<idx, index_type, vertex_id, Attributes...> const& G) {
    return G.parent();
  }

  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
            typename... Attributes, std::unsigned_integral VertexId>
  bool tag_invoke(is_local_index_tag,
    partitioned_index_adjacency_local_view<idx, index_type, vertex_id, Attributes...> const& G,
    VertexId v) {
    return G.is_local_index(static_cast<index_type>(v));
  }

  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
            typename... Attributes>
  partition_descriptor tag_invoke(partition_tag,
    partitioned_index_adjacency_local_view<idx, index_type, vertex_id, Attributes...> const& G) {
    return G.partition();
  }

  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
            typename... Attributes, std::unsigned_integral VertexId>
    requires std::convertible_to<VertexId, vertex_id>
  partition_descriptor tag_invoke(vertex_partition_tag,
    partitioned_index_adjacency_local_view<idx, index_type, vertex_id, Attributes...> const& G,
    VertexId v) {
    return vertex_partition(parent(G), v);
  }

  template <typename G>
  concept partitioned_graph_local_view =
    adjacency_list_graph<G> &&
    requires(G& g, std::size_t global_idx) {
      typename std::remove_reference_t<G>::graph_type;
      { parent(g) } -> partitioned_graph;
      { is_local_index(g, global_idx) } -> std::convertible_to<bool>;
    };

  template <typename G>
  concept partitioned_graph_with_local_view =
    partitioned_graph<G> &&
    requires(G& g, partition_t<G> partition) {
      { local_view(g, partition) } -> partitioned_graph_local_view;
    };


  ///**
  // * @brief Index biadjacency structure. This data structures stores bipartite graph in Compressed
  // * Sparse Row format. The underlying data structure is a structure of arrays for storage.
  // * Index_biadjacency is an biadjacency list represenation of a bipartite graph.
  // * A bipartite graph has two vertex partitions.
  // *
  // *
  // * @tparam idx The index type to indicate the type of the graph, can be either 0 or 1.
  // * @tparam index_type The data type used to represent a vertex index, required to be os unsigned
  // * integral type.
  // * @tparam vertex_id The data type used to represent a vertex ID, required to be os unsigned
  // * integral type.
  // * @tparam Attributes A variadic list of edge property types.
  // */
  //template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id,
  //          typename... Attributes>
  //class partitioned_index_biadjacency
  //  : public bipartite_graph_base,
  //    public hpx::partitioned_vector<
  //      index_type, indexed_struct_of_arrays<index_type, vertex_id, Attributes...>> {
  //  using base =
  //    hpx::partitioned_vector<index_type,
  //                            indexed_struct_of_arrays<index_type, vertex_id, Attributes...>>;

  //public:
  //  using value_type = index_type;
  //  using index_t = index_type;
  //  using vertex_id_type = vertex_id;
  //  using num_vertices_type = std::array<vertex_id_type, 2>;
  //  using num_edges_type = index_t;

  //  // The first index_t isn't considered an attribute.
  //  using attributes_t = std::tuple<Attributes...>;
  //  static constexpr std::size_t getNAttr() { return sizeof...(Attributes); }

  //  /**
  //   * @brief Constructor of index_biadjacency. Require the type of the graph to be bipartite.
  //   * Create an empty index_biadjacency.
  //   */
  //  partitioned_index_biadjacency(size_t N0 = 0, size_t N1 = 0, size_t M = 0)
  //    requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
  //    : bipartite_graph_base(N0, N1)
  //    , base(N0, M) {}
  //  /**
  //   * @brief Constructor of index_biadjacency. Require the type of the graph to be bipartite.
  //   * Create an empty index_biadjacency.
  //   */
  //  partitioned_index_biadjacency(std::array<size_t, 2> N, size_t M = 0)
  //    requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
  //    : bipartite_graph_base(N[idx], N[(idx + 1) % 2])
  //    , base(N[idx], M) {}

  //  template <class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
  //  partitioned_index_biadjacency(index_edge_list<vertex_id_type, bipartite_graph_base,
  //                                                directedness::directed, Attributes...>& A,
  //                                bool sort_biadjacency = false, ExecutionPolicy&& policy = {})
  //    requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
  //    : bipartite_graph_base(A.num_vertices()[idx], A.num_vertices()[(idx + 1) % 2])
  //    , base(A.num_vertices()[idx] + 1) {
  //    fill_biadjacency<idx>(A, *this, sort_biadjacency, policy);
  //  }

  //  template <class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
  //  partitioned_index_biadjacency(index_edge_list<vertex_id_type, bipartite_graph_base,
  //                                                directedness::undirected, Attributes...>& A,
  //                                bool sort_biadjacency = false, ExecutionPolicy&& policy = {})
  //    requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
  //    : bipartite_graph_base(A.num_vertices()[idx], A.num_vertices()[(idx + 1) % 2])
  //    , base(A.num_vertices()[idx] + 1) {
  //    fill_biadjacency<idx>(A, *this, sort_biadjacency, policy);
  //  }
  //  // customized move constructor
  //  partitioned_index_biadjacency(size_t N1, std::vector<vertex_id>&& indices,
  //                                std::vector<vertex_id>&& first_to_be,
  //                                std::vector<Attributes>&&... rest_to_be)
  //    requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
  //    : bipartite_graph_base(indices.size() - 1, N1)
  //    , base(std::move(indices), std::move(first_to_be), std::move(rest_to_be)...) {}
  //  partitioned_index_biadjacency(
  //    size_t N1, std::vector<vertex_id>&& indices,
  //    std::tuple<std::vector<vertex_id>, std::vector<Attributes>...>&& to_be_indexed)
  //    : bipartite_graph_base(indices.size() - 1, N1)
  //    , base(std::move(indices), std::move(to_be_indexed)) {}
  //  // customized copy constructor
  //  partitioned_index_biadjacency(size_t N1, const std::vector<vertex_id>& indices,
  //                                const std::vector<vertex_id>& first_to_be,
  //                                const std::vector<Attributes>&... rest_to_be)
  //    requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
  //    : bipartite_graph_base(indices.size() - 1, N1)
  //    , base(indices, first_to_be, rest_to_be...) {}
  //  partitioned_index_biadjacency(
  //    size_t N1, const std::vector<vertex_id>& indices,
  //    const std::tuple<std::vector<vertex_id>, std::vector<Attributes>...>& to_be_indexed)
  //    requires(std::is_same<bipartite_graph_base, bipartite_graph_base>::value)
  //    : bipartite_graph_base(indices.size() - 1, N1)
  //    , base(indices, to_be_indexed) {}

    //auto num_vertices() const { return vertex_cardinality; }
  //  num_edges_type num_edges() const { return base::to_be_indexed_.size(); };
  //  /**
  //   * @brief Serialize the partitioned_index_adjacency into binary file.
  //   *
  //   * @param outfile_name The output file name.
  //   */
  //  void serialize(const std::string& outfile_name) {
  //    std::ofstream out_file(outfile_name, std::ofstream::binary);
  //    bipartite_graph_base::serialize(out_file);
  //    base::serialize(out_file);
  //  }

  //  /**
  //   * @brief Deserialize the binary into partitioned_index_adjacency.
  //   *
  //   * @param infile_name The input file name.
  //   */
  //  void deserialize(const std::string& infile_name) {
  //    std::ifstream infile(infile_name, std::ifstream::binary);
  //    bipartite_graph_base::deserialize(infile);
  //    base::deserialize(infile);
  //  }
  //};

  //template <int idx, typename... Attributes>
  //using partitioned_biadjacency =
  //  partitioned_index_biadjacency<idx, default_index_t, default_vertex_id_type, Attributes...>;

  //template <int idx, edge_list_graph edge_list_t>
  //auto make_partitioned_biadjacency(edge_list_t& el) {
  //  return partitioned_biadjacency<idx>(el);
  //}


  //template <int idx, edge_list_c edge_list_t, std::unsigned_integral u_integral,
  //          class ExecutionPolicy = std::execution::parallel_unsequenced_policy>
  //auto make_partitioned_biadjacency(edge_list_t& el, u_integral n0, u_integral n1,
  //                                  directedness edge_directedness = directedness::directed,
  //                                  ExecutionPolicy&& policy = {}) {
  //  biadjacency<idx> adj(n0, n1);
  //  fill_biadjacency<idx>(el, adj, policy);
  //  return adj;
  //}


  //// template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
  //// typename... Attributes> auto num_vertices(const partitioned_index_adjacency<idx, index_type,
  //// vertex_id_type, Attributes...>& g) {
  ////   return g.num_vertices();
  //// }
  //// partitioned_index_adjacency_local_view num_vertices CPO
  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
            typename... Attributes>
  auto
  tag_invoke(const num_vertices_tag,
                  const partitioned_index_adjacency_local_view<idx, index_type, vertex_id_type,
                                                               Attributes...>& g) {
    return g.num_vertices();
  }
  // partitioned_index_adjacency_local_view degree CPO
  template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
            std::unsigned_integral lookup_type, typename... Attributes>
  auto
  tag_invoke(const degree_tag,
    const partitioned_index_adjacency_local_view<idx, index_type, vertex_id_type, Attributes...>& g,
             lookup_type i) {
    return g[i].size();
  }
  //// partitioned_index_adjacency degree CPO
  //template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
  //          typename... Attributes>
  //auto
  //tag_invoke(const degree_tag,
  //           const partitioned_index_adjacency<idx, index_type, vertex_id_type, Attributes...>& g,
  //           const typename partitioned_index_adjacency<idx, index_type, vertex_id_type,
  //                                                      Attributes...>::sub_view& v) {
  //  return v.size();
  //}
  //// partitioned_index_biadjacency num_vertices CPO
  //template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
  //          typename... Attributes>
  //auto
  //tag_invoke(const num_vertices_tag,
  //           const partitioned_index_biadjacency<idx, index_type, vertex_id_type, Attributes...>& g,
  //           int jdx = 0) {
  //  return g.num_vertices()[jdx];
  //}
  //// partitioned_index_biadjacency degree CPO
  //template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
  //          std::unsigned_integral lookup_type, typename... Attributes>
  //auto
  //tag_invoke(const degree_tag,
  //           const partitioned_index_biadjacency<idx, index_type, vertex_id_type, Attributes...>& g,
  //           lookup_type i) {
  //  return g[i].size();
  //}
  //// partitioned_index_biadjacency degree CPO
  //template <int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type,
  //          typename... Attributes>
  //auto
  //tag_invoke(const degree_tag,
  //           const partitioned_index_biadjacency<idx, index_type, vertex_id_type, Attributes...>& g,
  //           const typename partitioned_index_biadjacency<idx, index_type, vertex_id_type,
  //                                                        Attributes...>::sub_view& v) {
  //  return v.size();
  //}

  
  // New source/target CPOs
  template <typename Edge, int idx, std::unsigned_integral index_type,
            std::unsigned_integral vertex_id_type, typename... Attributes>
  auto tag_invoke(
    const source_tag,
    const partitioned_index_adjacency_local_view<idx, index_type, vertex_id_type, Attributes...>& g,
    Edge e) {
    return std::get<0>(e);
  }

  template <typename Edge, int idx, std::unsigned_integral index_type, std::unsigned_integral vertex_id_type, typename... Attributes>
  auto tag_invoke(const target_tag,
                  const partitioned_index_adjacency_local_view<idx, index_type, vertex_id_type,
                                                               Attributes...>& g, Edge e) {
    return std::get<0>(std::get<1>(e));
  }


// Sanity check, concepts should be satisfied
  template <typename G>
  concept has_member_begin = requires(G g) {
    { g.begin() } -> std::input_or_output_iterator;
  };

  template <typename G>
  concept has_member_end = requires(G g) {
    { g.end() } -> std::input_or_output_iterator;
  };

  static_assert(has_member_begin<partitioned_adjacency_local_view<0, int>>,
                "partitioned_adjacency_local_view does not have begin() member function");
    static_assert(has_member_end<partitioned_adjacency_local_view<0, int>>,
                "partitioned_adjacency_local_view does not have end() member function");
  static_assert(std::ranges::range<partitioned_adjacency_local_view<0, int>>,
                "partitioned_adjacency_local_view does not satisfy range concept");
    static_assert(adjacency_list_graph<partitioned_adjacency_local_view<0, int>>,
                "partitioned_adjacency_local_view does not satisfy adjacency_list_graph concept");

} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_ADJACENCY_HPP
