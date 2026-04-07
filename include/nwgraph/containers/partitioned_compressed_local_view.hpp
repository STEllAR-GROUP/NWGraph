/**
 * @file compressed.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * @authors
 *   Andrew Lumsdaine
 *   Luke D'Alessandro
 *   Kevin Deweese
 *   Krzysztof Drewniak
 *   Tony Liu
 *
 */

#ifndef NW_GRAPH_PARTITIONED_COMPRESSED_LOCAL_VIEW_HPP
#define NW_GRAPH_PARTITIONED_COMPRESSED_LOCAL_VIEW_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include "nwgraph/adaptors/splittable_range_adaptor.hpp"

#include "nwgraph/containers/partitioned_soa_local_view.hpp"
#include "nwgraph/util/partitioned_vector_local_partition_view.hpp"
#include "nwgraph/util/constant_iterator.hpp"

#include <hpx/algorithm.hpp>

namespace nw::graph {

  template <typename index_t, bool is_const = false, typename... Attributes>
  class partitioned_indexed_local_view_outer_iterator;

  template <typename index_t, typename... Attributes>
  class partitioned_indexed_struct_of_arrays_local_view {

    index_t first_;
    index_t last_;
    partition_descriptor partition_;


    nw::graph::util::partitioned_vector_local_partition_view<index_t> indices_;
    partitioned_struct_of_arrays_local_view<Attributes...> to_be_indexed_;

  public:
    template <typename T>
    using constant_iterator = nw::graph::util::constant_iterator<T>;
    // Modified to also store the source index
    using inner_iterator = hpx::util::zip_iterator<
      constant_iterator<index_t>,
      typename partitioned_struct_of_arrays_local_view<Attributes...>::iterator>;

    using const_inner_iterator = hpx::util::zip_iterator<
      constant_iterator<index_t>,
      typename partitioned_struct_of_arrays_local_view<Attributes...>::const_iterator>;

    using sub_view = nw::graph::splittable_range_adaptor<inner_iterator>;
    using const_sub_view = nw::graph::splittable_range_adaptor<const_inner_iterator>;

    static constexpr std::size_t getNAttr() { return sizeof...(Attributes); }

    partitioned_indexed_struct_of_arrays_local_view() = default;

    partitioned_indexed_struct_of_arrays_local_view(
      partitioned_indexed_struct_of_arrays<index_t, Attributes...>& parent_soa,
      partition_descriptor partition)
      : partition_(HPX_MOVE(partition))
      , indices_(parent_soa.indices_, partition_)
      , to_be_indexed_(
          parent_soa.to_be_indexed_,
          detail::aligned_partition(
            parent_soa.indices_,
            std::get<0>(
              static_cast<typename partitioned_struct_of_arrays<Attributes...>::base const&>(
                parent_soa.to_be_indexed_)),
            partition_)) {
      first_ = static_cast<index_t>(partition_.first_index());
      last_ = static_cast<index_t>(partition_.last_index());
      assert(indices_.is_local_index(first_));
      assert(last_ >= first_);
      assert(last_ <= indices_.last_index());
    }

    using const_outer_iterator =
      partitioned_indexed_local_view_outer_iterator<index_t, true, Attributes...>;
    using outer_iterator =
      partitioned_indexed_local_view_outer_iterator<index_t, false, Attributes...>;

    using iterator = partitioned_indexed_local_view_outer_iterator<index_t, false, Attributes...>;

    using value_type = iterator::value_type;
    using reference = iterator::reference;
    using size_type = std::size_t;
    using difference_type = iterator::difference_type;
    using pointer = iterator::pointer;

    using const_iterator = const_outer_iterator;
    using const_reference = const_iterator::reference;
    using const_pointer = const_iterator::pointer;

    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    // Remember, the local part of indices_ resides between first_ and last_ in global indexes
    iterator begin() { return {&indices_, &to_be_indexed_, first_, last_}; }
    const_iterator begin() const { return {&indices_, &to_be_indexed_, first_, last_}; }
    const_iterator cbegin() const { return {&indices_, &to_be_indexed_, first_, last_}; }
    iterator end() { return {&indices_, &to_be_indexed_, last_, last_}; }
    const_iterator end() const { return {&indices_, &to_be_indexed_, last_, last_}; }
    const_iterator cend() const { return {&indices_, &to_be_indexed_, last_, last_}; }

    /// Random access to the outer range (using global index).
    sub_view operator[](index_t i) { return *iterator{&indices_, &to_be_indexed_, i, last_}; }
    const_sub_view operator[](index_t i) const {
      return *iterator{&indices_, &to_be_indexed_, i, last_};
    }

    index_t size() const { return last_ - first_; }
    index_t max() const { return size() - 1; }

    auto& get_indices() { return indices_; }
    auto& get_to_be_indexed() { return to_be_indexed_; }

    auto const& get_indices() const { return indices_; }
    auto const& get_to_be_indexed() const { return to_be_indexed_; }

    
    bool is_local_index(index_t global_index) const {
      return indices_.is_local_index(global_index);
    }

    partition_descriptor partition() const { return partition_; }
  };

  template <typename index_t, bool is_const, typename... Attributes>
  class partitioned_indexed_local_view_outer_iterator {

  private:
    using isoa_t = partitioned_indexed_struct_of_arrays_local_view<index_t, Attributes...>;
    using constant_iterator = isoa_t::template constant_iterator<index_t>;
    using inner_iterator = isoa_t::inner_iterator;
    using sub_view = isoa_t::sub_view;
    using const_sub_view = isoa_t::const_sub_view;

    // Even though local views would be cheap to copy, they store a shared pointer, so let's avoid
    // storing copies of the views in the iterator itself, and instead store plain pointers.
    using indices_t =
      std::conditional_t<is_const, partitioned_vector_local_partition_view<index_t> const*,
                         partitioned_vector_local_partition_view<index_t>*>;
    using indexed_t =
      std::conditional_t<is_const, partitioned_struct_of_arrays_local_view<Attributes...> const*,
                         partitioned_struct_of_arrays_local_view<Attributes...>*>;

    indices_t indices_;
    indexed_t indexed_;
    index_t i_;
    index_t vertex_last_;

  public:
    using difference_type = std::make_signed_t<index_t>;
    using value_type = std::conditional_t<is_const, const_sub_view, sub_view>;
    using reference = value_type;
    using pointer = arrow_proxy<reference>;
    using iterator_category = std::random_access_iterator_tag;

    partitioned_indexed_local_view_outer_iterator() = default;

    partitioned_indexed_local_view_outer_iterator(
      indices_t indices, indexed_t indexed, index_t i, index_t vertex_last)
      : indices_(indices)
      , indexed_(indexed)
      , i_(i)
      , vertex_last_(vertex_last) {}

    partitioned_indexed_local_view_outer_iterator(
      partitioned_indexed_local_view_outer_iterator const&) = default;
    partitioned_indexed_local_view_outer_iterator(
      partitioned_indexed_local_view_outer_iterator<index_t, false, Attributes...> const& rhs)
      requires(is_const)
      : indices_(rhs.indices_)
      , indexed_(rhs.indexed_)
      , i_(rhs.i_)
      , vertex_last_(rhs.vertex_last_) {}

    partitioned_indexed_local_view_outer_iterator&
    operator=(partitioned_indexed_local_view_outer_iterator const&) = default;
    partitioned_indexed_local_view_outer_iterator& operator=(
      partitioned_indexed_local_view_outer_iterator<index_t, false, Attributes...> const& rhs)
      requires(is_const)
    {
      indices_ = rhs.indices_;
      indexed_ = rhs.indexed_;
      i_ = rhs.i_;
      vertex_last_ = rhs.vertex_last_;
      return *this;
    }

    partitioned_indexed_local_view_outer_iterator& operator++() {
      ++i_;
      return *this;
    }

    partitioned_indexed_local_view_outer_iterator operator++(int) {
      partitioned_indexed_local_view_outer_iterator tmp(*this);
      ++i_;
      return tmp;
      ;
    }

    partitioned_indexed_local_view_outer_iterator& operator--() {
      --i_;
      return *this;
    }

    partitioned_indexed_local_view_outer_iterator operator--(int) {
      partitioned_indexed_local_view_outer_iterator tmp(*this);
      --i_;
      return tmp;
    }

    partitioned_indexed_local_view_outer_iterator& operator+=(difference_type n) {
      i_ += n;
      return *this;
    }

    partitioned_indexed_local_view_outer_iterator& operator-=(difference_type n) {
      i_ -= n;
      return *this;
    }

    partitioned_indexed_local_view_outer_iterator operator+(difference_type n) const {
      return {indices_, indexed_, i_ + n, vertex_last_};
    }

    partitioned_indexed_local_view_outer_iterator operator-(difference_type n) const {
      return {indices_, indexed_, i_ - n, vertex_last_};
    }

    difference_type operator-(partitioned_indexed_local_view_outer_iterator const& b) const {
      return i_ - b.i_;
    }

    bool operator==(partitioned_indexed_local_view_outer_iterator const& b) const {
      return i_ == b.i_;
    }
    bool operator!=(partitioned_indexed_local_view_outer_iterator const& b) const {
      return i_ != b.i_;
    }
    bool operator<(partitioned_indexed_local_view_outer_iterator const& b) const {
      return i_ < b.i_;
    }
    bool operator>(partitioned_indexed_local_view_outer_iterator const& b) const {
      return i_ > b.i_;
    }
    bool operator<=(partitioned_indexed_local_view_outer_iterator const& b) const {
      return i_ <= b.i_;
    }
    bool operator>=(partitioned_indexed_local_view_outer_iterator const& b) const {
      return i_ >= b.i_;
    }

  private:
    auto indexed_current(index_t i) const { return indexed_->global_begin() + (*indices_)[i]; }

    // The edge range is constructed from {indexed_[i_], indexed_[i_+1]}
    // When i_+1 is local, we can directly use indexed_[i_+1]
    // But for the last local edge range, i_+1 may not be local, so in that
    // case the last edge range is terminated by indexed_.end(), essentially
    // containing all remaining local edges.
    auto indexed_next(index_t i) const {
      // If i+1 is the last local index, return the end of indexed_
      if (i + 1 == vertex_last_) {
        return indexed_->end();
      }
      else {
        return indexed_->global_begin() + (*indices_)[i + 1];
      }
    }

  public:
    reference operator*() {
      assert(indices_->is_local_index(i_));
      constant_iterator c_it{i_};
      inner_iterator start{c_it, indexed_current(i_)};
      inner_iterator end{c_it, indexed_next(i_)};
      return {start, end};
    }
    reference operator*() const {
      assert(indices_->is_local_index(i_));
      constant_iterator c_it{i_};
      inner_iterator start{c_it, indexed_current(i_)};
      inner_iterator end{c_it, indexed_next(i_)};
      return {start, end};
    }

    pointer operator->() { return {**this}; }
    pointer operator->() const { return {**this}; }

    reference operator[](index_t n) {
      assert(indices_->is_local_index(i_ + n));
      constant_iterator c_it{i_+n};
      inner_iterator start{c_it, indexed_current(i_ + n)};
      inner_iterator end{c_it, indexed_next(i_ + n)};
      return {start, end};
    }
    reference operator[](index_t n) const {
      assert(indices_->is_local_index(i_ + n));
      constant_iterator c_it{i_ + n};
      inner_iterator start{c_it, indexed_current(i_ + n)};
      inner_iterator end{c_it, indexed_next(i_ + n)};
      return {start, end};
    }

    // auto index() { return indices_.iter_from_global_index(i_); }
    // auto index() const { return indices_.iter_from_global_index(i_); }

    index_t index() const { return i_; }; // Maybe this is more sensible
  };

  // template <typename index_t, typename... Attributes>
  // auto operator+(std::iter_difference_t<typename partitioned_indexed_struct_of_arrays_local_view<
  //                  index_t, Attributes...>::outer_iterator>
  //                  n,
  //                typename partitioned_indexed_struct_of_arrays_local_view<
  //                  index_t, Attributes...>::outer_iterator const& i) {
  //   return i + n;
  // }

} // namespace nw::graph

#endif // NW_GRAPH_PARTITIONED_COMPRESSED_LOCAL_VIEW_HPP
