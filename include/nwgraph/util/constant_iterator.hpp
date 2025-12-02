#ifndef NW_GRAPH_UTIL_CONSTANT_ITERATOR_HPP
#define NW_GRAPH_UTIL_CONSTANT_ITERATOR_HPP

namespace nw::graph::util {

  // An iterator that generates a constant value
  template <typename T>
  class constant_iterator {

  public:
    constant_iterator() = default;

    explicit constant_iterator(T value)
      : value_(value) {}

    constant_iterator& operator++() { return *this; }
    constant_iterator operator++(int) { return *this; }
    constant_iterator& operator--() { return *this; }
    constant_iterator operator--(int) { return *this; }

    constant_iterator operator+(constant_iterator const&) const { return *this; }
    std::ptrdiff_t operator-(constant_iterator const&) const { return 0; }

    constant_iterator& operator+=(std::ptrdiff_t) { return *this; }
    constant_iterator& operator-=(std::ptrdiff_t) { return *this; }
    constant_iterator operator+(std::ptrdiff_t) const { return *this; }
    constant_iterator operator-(std::ptrdiff_t) const { return *this; }

    // Never compare equal
    // This is to allow zip_iterator to terminate based on other iterators
    bool operator==(constant_iterator const&) const { return false; }
    auto operator<=>(constant_iterator const&) const = default;

    T operator*() const { return value_; }
    T operator[](std::ptrdiff_t) const { return value_; }


    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using reference = T;
    using iterator_category = std::random_access_iterator_tag;

  private:
    T value_;
  };
} // namespace nw::graph::util

#endif // NW_GRAPH_UTIL_CONSTANT_ITERATOR_HPP