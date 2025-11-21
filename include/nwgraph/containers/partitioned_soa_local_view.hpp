/**
 * @file partitioned_soa_local_view.hpp
 */

#ifndef NW_GRAPH_PARTITIONED_SOA_LOCAL_VIEW_HPP
#define NW_GRAPH_PARTITIONED_SOA_LOCAL_VIEW_HPP

#ifndef NWGRAPH_HAVE_HPX
#error "This file requires using HPX as a backend for NWGraph"
#endif

#include <cassert>

#include <fstream>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <istream>
#include <iterator>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "nwgraph/containers/soa.hpp"
#include "nwgraph/util/partitioned_vector_local_partition_view.hpp"

#include <algorithm>
#include <execution>

#include <hpx/include/partitioned_vector.hpp>
#include <hpx/include/runtime.hpp>
#include <hpx/include/serialization.hpp>

namespace nw::graph {

  using nw::graph::util::partitioned_vector_local_partition_view;

  template <class... Attributes>
  struct partitioned_struct_of_arrays_local_view
    : std::tuple<partitioned_vector_local_partition_view<Attributes>...> {
    using base = std::tuple<partitioned_vector_local_partition_view<Attributes>...>;

    std::size_t loc_begin_;
    std::size_t loc_end_;

  public:
    template <bool is_const = false>
    class soa_iterator {
      friend class soa_iterator<!is_const>;

      using soa_t = std::conditional_t<is_const, partitioned_struct_of_arrays_local_view const,
                                       partitioned_struct_of_arrays_local_view>;
      std::size_t i_{0};
      soa_t* soa_{nullptr}; // Can this be a reference?

    public:
      using _val_type =
        std::tuple<typename partitioned_vector_local_partition_view<Attributes>::value_type...>;
      using _cval_type = std::tuple<
        typename partitioned_vector_local_partition_view<Attributes>::value_type const...>;

      using value_type = std::conditional_t<is_const, _cval_type, _val_type>;

      using _ref_type =
        std::tuple<typename partitioned_vector_local_partition_view<Attributes>::reference...>;
      using _cref_type = std::tuple<
        typename partitioned_vector_local_partition_view<Attributes>::const_reference...>;

      using reference = std::conditional_t<is_const, _cref_type, _ref_type>;

      using difference_type = std::ptrdiff_t;

      using pointer = arrow_proxy<reference>; // How is this used?
      using iterator_category = std::random_access_iterator_tag;

      soa_iterator() = default;

      explicit soa_iterator(soa_t* soa, std::size_t i = 0)
        : i_(i)
        , soa_(soa) {}

      soa_iterator(soa_iterator const&) = default;
      explicit soa_iterator(soa_iterator<false> const& b)
        requires(is_const)
        : i_(b.i_)
        , soa_(b.soa_) {}

      soa_iterator& operator=(soa_iterator const&) = default;

      soa_iterator& operator=(soa_iterator<false> const& b)
        requires(is_const)
      {
        i_ = b.i_;
        soa_ = b.soa_;
        return *this;
      }

      // This also compares the soa_ pointers, which I think is only necessary
      // for equality comparison, but not for ordering. We don't mind unless
      // it's performance critical.
      bool operator==(soa_iterator const&) const = default;
      auto operator<=>(soa_iterator const&) const = default;

      soa_iterator operator++(int) { return soa_iterator{soa_, i_++}; }

      soa_iterator operator--(int) { return soa_iterator{soa_, i_--}; }

      soa_iterator& operator++() {
        ++i_;
        return *this;
      }

      soa_iterator& operator--() {
        --i_;
        return *this;
      }

      soa_iterator& operator+=(std::ptrdiff_t n) {
        i_ += n;
        return *this;
      }

      soa_iterator& operator-=(std::ptrdiff_t n) {
        i_ -= n;
        return *this;
      }

      soa_iterator operator+(std::ptrdiff_t n) const { return soa_iterator(soa_, i_ + n); }

      soa_iterator operator-(std::ptrdiff_t n) const { return soa_iterator(soa_, i_ - n); }

      std::ptrdiff_t operator-(soa_iterator const& b) const { return i_ - b.i_; }

      friend soa_iterator operator+(std::ptrdiff_t n, soa_iterator i) { return i + n; }

      friend soa_iterator operator-(std::ptrdiff_t n, soa_iterator i) { return i - n; }

      decltype(auto) operator*() const {
        if constexpr (is_const) {
          return std::apply([this]<class... Vectors>(Vectors&&... v)
                            { return reference(std::forward<Vectors>(v)[i_]...); },
                            const_cast<base&>(static_cast<base const&>(*soa_)));
        }
        else {
          return std::apply([this]<class... Vectors>(Vectors&&... v)
                            { return reference(std::forward<Vectors>(v)[i_]...); },
                            static_cast<base&>(*soa_));
        }
      }

      decltype(auto) operator[](std::ptrdiff_t n) const {
        if constexpr (is_const) {
          return std::apply([this, n]<class... Vectors>(Vectors&&... v)
                            { return reference(std::forward<Vectors>(v)[i_ + n]...); },
                            const_cast<base&>(static_cast<base const&>(*soa_)));
        }
        else {
          return std::apply([this, n]<class... Vectors>(Vectors&&... v)
                            { return reference(std::forward<Vectors>(v)[i_ + n]...); },
                            static_cast<base&>(*soa_));
        }
      }

      // How is this used?
      pointer operator->() const { return {**this}; }

      pointer operator->() { return {**this}; }

      auto index() const { return i_; }
    };

    using iterator = soa_iterator<false>;

    using value_type = typename iterator::value_type;
    using reference = typename iterator::reference;
    using size_type = std::size_t;
    using difference_type = typename iterator::difference_type;
    using pointer = typename iterator::pointer;

    using const_iterator = soa_iterator<true>;
    using const_reference = typename const_iterator::reference;
    using const_pointer = typename const_iterator::pointer;

    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    partitioned_struct_of_arrays_local_view() = default;

  private:
    // Helper to go from tuple<pv> to tuple<pv_local_view>
    base _constructor_tuple_helper(partitioned_struct_of_arrays<Attributes...>& soa,
                                   std::size_t partnum) {
      return std::apply(
        [&](auto&... pvec)
        {
          return std::make_tuple(
            partitioned_vector_local_partition_view(pvec, partnum)...);
        },
        static_cast<typename partitioned_struct_of_arrays<Attributes...>::base&>(soa));
    }

  public:
    explicit partitioned_struct_of_arrays_local_view(
      partitioned_struct_of_arrays<Attributes...>& soa, std::size_t partnum)
      : base(_constructor_tuple_helper(soa, partnum)) {}

    iterator begin() {
      std::size_t i = std::get<0>(static_cast<base&>(*this)).first_index();
      return iterator(this, i);
    }

    // Needed for begin()+global_idx access pattern
    iterator global_begin() { return iterator(this, 0); }

    const_iterator begin() const {
      std::size_t i = std::get<0>(static_cast<base const&>(*this)).first_index();
      return const_iterator(this, i);
    }

    const_iterator cbegin() const {
      std::size_t i = std::get<0>(static_cast<base const&>(*this)).first_index();
      return const_iterator(this, i);
    }

    iterator end() { return begin() + size(); }

    const_iterator end() const { return begin() + size(); }

    const_iterator cend() const { return begin() + size(); }

    reference operator[](std::size_t i) {
      return std::apply([&](auto&... r)
                        { return std::forward_as_tuple(std::forward<decltype(r)>(r)[i]...); },
                        static_cast<base&>(*this));
    }

    const_reference operator[](std::size_t i) const {
      return std::apply([&](auto&... r)
                        { return std::forward_as_tuple(std::forward<decltype(r)>(r)[i]...); },
                        static_cast<base const&>(*this));
    }

    constexpr pointer data() noexcept {
      return std::apply([&](auto&... p)
                        { return std::forward_as_tuple(std::forward<decltype(p)>(p).data()...); },
                        static_cast<base&>(*this));
    }

    constexpr const_pointer data() const noexcept {
      return std::apply([&](auto&... p)
                        { return std::forward_as_tuple(std::forward<decltype(p)>(p).data()...); },
                        static_cast<base const&>(*this));
    }

    [[nodiscard]] size_t size() const {
      return std::get<0>(static_cast<base const&>(*this)).size();
    }

    bool is_local_index(std::size_t global_idx) const {
      return std::get<0>(static_cast<base const&>(*this)).is_local_index(global_idx);
    }
  };
} // namespace nw::graph


// Do we need this?
namespace std {
  template <class... Attributes>
  class tuple_size<nw::graph::partitioned_struct_of_arrays_local_view<Attributes...>>
    : public std::integral_constant<std::size_t, sizeof...(Attributes)> {};

} // namespace std

/// NB: technically we're supposed to be using `iter_swap` here on the
/// struct_of_array iterator type, but I can't figure out how to do this.

#include "nwgraph/util/tuple_hack.hpp"

#endif // NW_GRAPH_PARTITIONED_SOA_LOCAL_VIEW_HPP
