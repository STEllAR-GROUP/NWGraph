
#ifndef NW_GRAPH_PARTITIONED_VECTOR_LOCAL_PARTITION_VIEW
#define NW_GRAPH_PARTITIONED_VECTOR_LOCAL_PARTITION_VIEW

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/partitioned_adjacency.hpp"

#include <hpx/components/get_ptr.hpp>
#include <hpx/include/partitioned_vector.hpp>
#include <algorithm>
#include <memory>
#include <vector>

namespace nw::graph::util {

  // partitioned_vector_local_partition_view
  // It's purpose is to allow direct, cheap access to the data of a single, local partition
  // using the global indexes. It should also be able to answer whether an index belongs
  // to it quite quickly.

  template <typename T>
  class partitioned_vector_local_partition_view {

    friend class hpx::partitioned_vector<T>;

    // Useful to keep a (non-owning) referece to the global (parent) partitioned vector
    hpx::partitioned_vector<T>* parent_ = nullptr;

    // I am a stubborn person
    using partitioned_vector_server = decltype(std::declval<hpx::partitioned_vector<T>>()
                                                 .partitions()[0]
                                                 .local_data_)::element_type;

    // Could alternatively store direct reference to internal Data;
    std::shared_ptr<partitioned_vector_server> data_;


    // Global index of the first element in this partition
    std::size_t first_;
    std::size_t size_;
    partition_descriptor partition_;

    // friend class hpx::serialization::access;

    // void serialize(hpx::serialization::input_archive& ar, unsigned) {
    //   ar >> parent_ >> first_ >> size_ >> data_;
    // }

    // void serialize(hpx::serialization::output_archive& ar, unsigned) const {
    //   ar << parent_ << first_ << size_ << data_;
    // }

    std::size_t local_index_from_global(std::size_t global_idx) const {
      HPX_ASSERT(is_local_index(global_idx));
      return global_idx - first_;
    }

  public:
    using value_type = T;
    using reference = T&;
    using const_reference = T const&;

    // partitioned_vector_local_partition_view() = default;

    partitioned_vector_local_partition_view(
      hpx::partitioned_vector<T>& pv, partition_descriptor partition)
      : parent_(&pv)
      , partition_(HPX_MOVE(partition)) {
      auto component_id = detail::partition_id(partition_);
      if (detail::find_partition_by_id(pv, partition_) == pv.segment_end()) {
        component_id = detail::find_matching_partition(pv, partition_)->get_id();
      }
      first_ = partition_.first_index();
      size_ = partition_.size();
      data_ = hpx::get_ptr<partitioned_vector_server>(hpx::launch::sync, component_id);
      HPX_ASSERT(data_);
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
    // TODO: maybe return Data::iterator instead?
    T* begin() { return data_->get_data().data(); }
    T const* begin() const { return data_->get_data().data(); }
    T* end() { return data_->get_data().data() + size_; }
    T const* end() const { return data_->get_data().data() + size_; }

    T* iter_from_global_index(std::size_t global_idx) {
      std::size_t local_idx = local_index_from_global(global_idx);
      return data_->get_data().data() + local_idx;
    }

    std::size_t size() const { return size_; }
    std::size_t first_index() const { return first_; }
    std::size_t last_index() const { return first_ + size_; }
    partition_descriptor partition() const { return partition_; }

    hpx::partitioned_vector<T>& parent() {
        assert(parent_);
        return *parent_; 
    }

    hpx::partitioned_vector<T> const& parent() const {
        assert(parent_);
        return *parent_;
    }
  };

} // namespace nw::graph::util

#endif // ! NW_GRAPH_PARTITIONED_VECTOR_LOCAL_PARTITION_VIEW
