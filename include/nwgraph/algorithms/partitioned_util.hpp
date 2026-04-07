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
    out_degree_count_packet(
      hpx::partitioned_vector<vertex_id_t<std::remove_reference_t<Graph>>> degrees,
      std::vector<vertex_id_t<std::remove_reference_t<Graph>>>&& incoming_packet) {

      for (auto&& v_dest : incoming_packet) {
        auto deg_iter = degrees.get_local_iterator(v_dest).local();
        (*deg_iter)++;
      }
    }

    template <typename Graph>
    struct out_degree_count_action
      : hpx::actions::action<decltype(&out_degree_count_packet<Graph>),
                             &out_degree_count_packet<Graph>, out_degree_count_action<Graph>> {};


    struct out_degree_count : hpx::parallel::detail::algorithm<out_degree_count, int> {

      constexpr out_degree_count() noexcept
        : hpx::parallel::detail::algorithm<out_degree_count, int>("out_degree_count") {}


      template <typename ExPolicy, typename Graph>
      static int
      sequential(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index,
                             hpx::partitioned_vector<vertex_id_t<std::remove_reference_t<Graph>>> degrees) {

        // Assume adjacency list contains incoming edges

        hpx::id_type this_locality_id = hpx::find_here();

        auto first = G.begin() + first_index;
        auto last = G.begin() + last_index;

        auto deg_iter = degrees.get_local_iterator(first_index).local();

        using vertex_id_type = vertex_id_t<std::remove_reference_t<Graph>>;

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

        return 0;
      }

      template <typename ExPolicy, typename Graph>
      static int
      parallel(ExPolicy&& policy, Graph G, const size_t first_index,
                             const size_t last_index,
                             hpx::partitioned_vector<vertex_id_t<std::remove_reference_t<Graph>>> degrees) {
      // boop
        return 0;
      }
    };



    

    struct in_degree_count : hpx::parallel::detail::algorithm<in_degree_count, int> {
      static constexpr bool partition_aware = true;

      constexpr in_degree_count() noexcept
        : hpx::parallel::detail::algorithm<in_degree_count, int>("in_degree_count") {}


      template <typename ExPolicy, typename Graph>
      static int sequential(ExPolicy&& policy, Graph G, partition_descriptor partition,
                             hpx::partitioned_vector<vertex_id_t<std::remove_reference_t<Graph>>> degrees) {

        auto G_loc = local_view(G, partition);
        auto degrees_loc = local_view(degrees, partition);

        auto u = static_cast<vertex_id_t<std::remove_reference_t<Graph>>>(partition.first_index());
        for (auto v_it = G_loc.begin(); v_it != G_loc.end(); ++v_it, ++u) {
          degrees_loc[u] = v_it->size();
        }

        return 0;

      }

      template <typename ExPolicy, typename Graph>
      static int parallel(ExPolicy&& policy, Graph G, const size_t first_index,
                           const size_t last_index,
                           hpx::partitioned_vector<vertex_id_t<std::remove_reference_t<Graph>>> degrees) {
        // boop
          return 0;
      }
    };


  } // namespace detail


  template <partitioned_algorithm_graph Graph>
  auto partitioned_degrees(Graph& G) {

    using vertex_id_type = vertex_id_t<std::remove_reference_t<Graph>>;

    auto degrees = make_copartitioned_vector<vertex_id_type>(G);

    degrees.register_as("degrees"); // TODO: Do we need a unique name on each invocation?
    

    //partitioned_algorithm<detail::out_degree_count>(hpx::execution::seq, G, hpx::ref(degrees));
    auto futures = partitioned_algorithm<detail::in_degree_count>(hpx::execution::seq, G, hpx::ref(degrees));
    hpx::wait_all(futures);

    return degrees;
  }











    namespace detail {
    struct avg_degree_per_partition
      : hpx::parallel::detail::algorithm<avg_degree_per_partition, double> {
      static constexpr bool partition_aware = true;

      constexpr avg_degree_per_partition() noexcept
        : hpx::parallel::detail::algorithm<avg_degree_per_partition, double>(
            "avg_degree_per_partition") {}


      template <typename ExPolicy, typename Graph>
      static double sequential(ExPolicy&& policy, Graph G, partition_descriptor partition) {

        auto G_loc = local_view(G, partition);

        double count = 0;
        double part_size = G_loc.size();

        // for each v in G do
        for (auto v_it = G_loc.begin(); v_it != G_loc.end(); ++v_it) {
          count += v_it->size();
        }

        return count / part_size;
      }

      template <typename ExPolicy, typename Graph>
      static float parallel(ExPolicy&& policy, Graph G, const size_t first_index,
                            const size_t last_index) {
        // boop
        return 0;
      }
    };
  } // namespace detail


    template <partitioned_algorithm_graph Graph>
    std::vector<double> partitioned_avg_degree_per_partition(Graph& G) {
      
        auto results = partitioned_algorithm<detail::avg_degree_per_partition>(hpx::execution::seq, G);
      
      // unpack futures
      std::vector<double> degrees;
      for (auto&& f : results) {
        degrees.push_back(f.get());
      }
      
      return degrees;
    }


    namespace detail {
      struct avg_remote_degree_per_partition
        : hpx::parallel::detail::algorithm<avg_remote_degree_per_partition, double> {
        static constexpr bool partition_aware = true;

        constexpr avg_remote_degree_per_partition() noexcept
          : hpx::parallel::detail::algorithm<avg_remote_degree_per_partition, double>(
              "avg_remote_degree_per_partition") {}

        template <typename ExPolicy, typename Graph>
        static double sequential(ExPolicy&& policy, Graph G, partition_descriptor partition) {
          hpx::id_type this_locality_id = hpx::find_here();
          auto G_loc = local_view(G, partition);
          double count = 0;
          double part_size = G_loc.size();
          // for each v in G do
          for (auto v_it = G_loc.begin(); v_it != G_loc.end(); ++v_it) {
            for (auto&& neighbor : *v_it) {
              auto v = target(G_loc, neighbor);
              if (!is_local(partition, v)) {
                count++;
              }
            }
          }
          return count / part_size;
        }

        template <typename ExPolicy, typename Graph>
        static float parallel(ExPolicy&& policy, Graph G, const size_t first_index,
                              const size_t last_index) {
          // boop
          return 0;
        }
      };
    } // namespace detail

    template <partitioned_algorithm_graph Graph>
    std::vector<double> partitioned_avg_remote_degree_per_partition(Graph& G) {
      auto results = partitioned_algorithm<detail::avg_remote_degree_per_partition>(hpx::execution::seq, G);
      // unpack futures
      std::vector<double> degrees;
      for (auto&& f : results) {
        degrees.push_back(f.get());
      }
      return degrees;
    }

} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_UTIL_HPP
