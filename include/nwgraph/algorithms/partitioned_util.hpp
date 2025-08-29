#ifndef NW_GRAPH_PARTITIONED_UTIL_HPP
#define NW_GRAPH_PARTITIONED_UTIL_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/algorithms/partitioned_algorithm.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <vector>

#include <hpx/async_combinators/wait_all.hpp>
#include <hpx/executors/execution_policy.hpp>
#include <hpx/include/partitioned_vector.hpp>
#include <hpx/include/partitioned_vector_predef.hpp>
#include <hpx/parallel/segmented_algorithms/detail/dispatch.hpp>
#include <hpx/parallel/util/detail/algorithm_result.hpp>


namespace nw::graph {

  namespace detail {

    template <typename Graph>
    static void
    out_degree_count_packet(hpx::partitioned_vector<typename Graph::vertex_id_type> degrees,
                            std::vector<typename Graph::vertex_id_type>&& incoming_packet) {

      for (auto&& v_dest : incoming_packet) {
        auto deg_iter = degrees.get_local_iterator(v_dest).local();
        (*deg_iter)++;
      }
    }

    template <typename Graph>
    struct out_degree_count_action
      : hpx::actions::action<decltype(&out_degree_count_packet<Graph>),
                             &out_degree_count_packet<Graph>, out_degree_count_action<Graph>> {};


    struct out_degree_count : hpx::parallel::detail::algorithm<out_degree_count> {

      constexpr out_degree_count() noexcept
        : hpx::parallel::detail::algorithm<out_degree_count>("out_degree_count") {}


      template <typename ExPolicy, typename Graph>
      static void sequential(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index,
                             hpx::partitioned_vector<typename Graph::vertex_id_type> degrees) {

        // Assume adjacency list contains incoming edges

        hpx::id_type this_locality_id = hpx::find_here();

        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;

        auto deg_iter = degrees.get_local_iterator(first_index).local();

        using vertex_id_type = typename Graph::vertex_id_type;

        std::map<hpx::id_type, std::vector<vertex_id_type>> outgoing_packets;

        // for each v in G do
        for (auto v_it = first; v_it != last; ++v_it) {

          auto neighbor_range = *v_it;
          for (auto neighbor : neighbor_range) {

            vertex_id_type v = target(G, neighbor);

            auto v_loc_id = vertex_locality(G, v);
            if (v_loc_id == this_locality_id) {
              // send to elt's locality
              outgoing_packets[v_loc_id].push_back(v);
            }
            else {
              // handle things locally
              auto deg_iter = degrees.get_local_iterator(v).local();
              (*deg_iter)++;
            }
          }
        }

        // Now send counts to other localities
        std::vector<hpx::future<void>> remote_ops;
        remote_ops.reserve(outgoing_packets.size());

        for (auto&& [id, packet] : outgoing_packets) {
          using act_t = out_degree_count_action<Graph>;
          remote_ops.push_back(hpx::async<act_t>(id, hpx::ref(degrees), std::move(packet)));
        }

        hpx::wait_all(remote_ops);
      }

      template <typename ExPolicy, typename Graph>
      static void parallel(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index,
                             hpx::partitioned_vector<typename Graph::vertex_id_type> degrees) {
      // boop
      }
    };



    

    struct in_degree_count : hpx::parallel::detail::algorithm<in_degree_count> {

      constexpr in_degree_count() noexcept
        : hpx::parallel::detail::algorithm<in_degree_count>("in_degree_count") {}


      template <typename ExPolicy, typename Graph>
      static void sequential(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index,
                             hpx::partitioned_vector<typename Graph::vertex_id_type> degrees) {

        // Assume adjacency list contains incoming edges

        hpx::id_type this_locality_id = hpx::find_here();

        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;

        auto deg_iter = degrees.get_local_iterator(first_index).local();

        using vertex_id_type = typename Graph::vertex_id_type;

        // for each v in G do
        for (auto v_it = first; v_it != last; ++v_it, ++deg_iter) {
          *deg_iter = v_it->size();
        }

      }

      template <typename ExPolicy, typename Graph>
      static void parallel(ExPolicy&& policy, Graph G, const size_t first_index,
                           const size_t last_index,
                           hpx::partitioned_vector<typename Graph::vertex_id_type> degrees) {
        // boop
      }
    };


  } // namespace detail


  auto partitioned_degrees(auto& G) {

    using vertex_id_t = typename std::decay_t<decltype(G)>::vertex_id_type;

    hpx::partitioned_vector<vertex_id_t> degrees(
      G.size(),
      hpx::explicit_container_layout(G.indices_.get_partition_sizes(),
                                     G.indices_.get_partition_localities()));

    degrees.register_as("degrees"); // TODO: Do we need a unique name on each invocation?
    

    //partitioned_algorithm<detail::out_degree_count>(hpx::execution::seq, G, hpx::ref(degrees));
    auto futures = partitioned_algorithm<detail::in_degree_count>(hpx::execution::seq, G, hpx::ref(degrees));
    hpx::wait_all(futures);

    return degrees;
  }

} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_UTIL_HPP
