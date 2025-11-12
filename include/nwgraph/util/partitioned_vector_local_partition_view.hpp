
#ifndef NW_GRAPH_PARTITIONED_VECTOR_LOCAL_PARTITION_VIEW
#define NW_GRAPH_PARTITIONED_VECTOR_LOCAL_PARTITION_VIEW

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include <hpx/include/partitioned_vector_predef.hpp>
#include <memory>


namespace nw::graph::util {

  // partitioned_vector_local_partition_view
  // It's purpose is to allow direct, cheap access to the data of a single, local partition
  // using the global indexes. It should also be able to answer whether an index belongs
  // to it quite quickly.

  template <typename T, typename Data>
  class partitioned_vector_local_partition_view {

    using partitioned_vector_server = hpx::server::partitioned_vector<T, Data>;

    std::shared_ptr<partitioned_vector_server> data_;

    // or, could store direct reference to internal Data
    // Data& data_;

    // Global index of the first element in this partition
    std::size_t first_;
    std::size_t size_;

    std::size_t local_index_from_global(std::size_t global_idx) const {
      HPX_ASSERT(is_local_index(global_idx));
      return global_idx - first_;
    }

  public:
    partitioned_vector_local_partition_view(hpx::partitioned_vector<T, Data> const& pv,
                                            std::size_t partnum) {
      HPX_ASSERT(partnum < pv.partitions().size());
      auto& partition = pv.partitions()[partnum];
      first_ = partition.first_;
      size_ = partition.size_;
      data_ = partition.local_data_;
    }

    bool is_local_index(std::size_t global_idx) const {
      return global_idx >= first_ && global_idx < first_ + size_;
    }

    /* [] operator */

    T& operator[](std::size_t global_idx) {
      std::size_t local_idx = local_index_from_global(global_idx);
      return data_->get_data()[local_idx];
      // or
      // return data_[local_idx];
    }

    T const& operator[](std::size_t global_idx) const {
      std::size_t local_idx = local_index_from_global(global_idx);
      return data_->get_data()[local_idx];
      // or
      // return data_[local_idx];
    }


    /* Iteration methods*/
    T* begin() { return data_->get_data().data(); }
    T const* begin() const { return data_->get_data().data(); }
    T* end() { return data_->get_data().data() + size_; }
    T const* end() const { return data_->get_data().data() + size_; }

    std::size_t size() const { return size_; }
  };

} // namespace nw::graph::util

#endif // ! NW_GRAPH_PARTITIONED_VECTOR_LOCAL_PARTITION_VIEW
