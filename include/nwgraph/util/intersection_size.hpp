/**
 * @file intersection_size.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * @authors
 *   Andrew Lumsdaine
 *   Luke D'Alessandro
 *
 */

#ifndef NW_GRAPH_INTERSECTION_SIZE_HPP
#define NW_GRAPH_INTERSECTION_SIZE_HPP

#if defined(CL_SYCL_LANGUAGE_VERSION)
#include <dpstd/algorithm>
#include <dpstd/execution>
#include <dpstd/numeric>
#elif defined(NWGRAPH_HAVE_HPX)
#include <hpx/algorithm.hpp>
#include <hpx/execution.hpp>
#include <hpx/numeric.hpp>
#else
#include <algorithm>
#include <execution>
#include <numeric>
#endif
#include <type_traits>

namespace nw {
  namespace graph {
    namespace detail {
      template <class It1, class It2, class Compare>
      std::size_t intersection_size_seq(It1 i, It1&& ie, It2 j, It2&& je, Compare cmp) {
        std::size_t n = 0;
        while (i != ie && j != je) {
          if (cmp(*i, *j)) {
            ++i;
          }
          else if (cmp(*j, *i)) {
            ++j;
          }
          else {
            ++n;
            ++i;
            ++j;
          }
        }
        return n;
      }

      template <class ExecutionPolicy, class It1, class It2, class Compare>
      std::size_t intersection_size_impl(ExecutionPolicy&& ep, It1 i, It1&& ie, It2 j, It2&& je,
                                         Compare cmp) {
        // Use our own trivial loop for the intersection size when the execution
        // policy is sequential, otherwise rely on std::set_intersection.
        //
        // @todo We really don't need set intersection. You'd hope that it would be
        //       efficient with the output counter, but it just isn't. Parallelizing
        //       the intersection size seems non-trivial though.
        if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
                                     std::execution::sequenced_policy>) {
          return detail::intersection_size_seq(i, std::forward<It1>(ie), j, std::forward<It2>(je),
                                               cmp);
        }
        else {
          return std::size_t(std::set_intersection(
            std::forward<ExecutionPolicy>(ep), std::forward<It1>(i), std::forward<It1>(ie),
            std::forward<It2>(j), std::forward<It2>(je), nw::graph::counter{}, cmp));
        }
      }

    } // namespace detail


    /// Basic helper used for all of the inner set intersections.
    ///
    /// This wraps `std::set_intersection` to produce the size of the set rather
    /// than the set itself, and also handles the fact that our iterator value types
    /// are tuples where we only care about the first element for ordering.
    ///
    /// @tparam ExecutionPolicy The type of the parallel execution policy.
    /// @tparam           It1 The type of the first iterator.
    /// @tparam           It2 The type of the second iterator.
    ///
    /// @param           ep The parallel execution policy.
    /// @param            i The beginning of the first range.
    /// @param           ie The end of the first range.
    /// @param            j The beginning of the second range.
    /// @param           je The end of the second range.
    ///
    /// @returns            The size of the intersected set.
    template <class ExecutionPolicy, std::forward_iterator It1, std::forward_iterator It2>
    std::size_t intersection_size(ExecutionPolicy&& ep, It1 i, It1&& ie, It2 j, It2&& je) {
      static constexpr auto lt = [](auto&& x, auto&& y) { return std::get<0>(x) < std::get<0>(y); };
      return detail::intersection_size_impl(std::forward<ExecutionPolicy>(ep), std::forward<It1>(i),
                                            std::forward<It1>(ie), std::forward<It2>(j),
                                            std::forward<It2>(je), lt);
    }

    // With custom comp function
    template <class ExecutionPolicy, std::forward_iterator It1, std::forward_iterator It2,
              class Compare>
    std::size_t intersection_size(ExecutionPolicy&& ep, It1 i, It1&& ie, It2 j, It2&& je,
                                  Compare cmp) {
      return detail::intersection_size_impl(std::forward<ExecutionPolicy>(ep), std::forward<It1>(i),
                                            std::forward<It1>(ie), std::forward<It2>(j),
                                            std::forward<It2>(je), cmp);
    }

    /// A convenience overload for `intersection_size`.
    ///
    /// This overload takes two ranges and an execution policy, and forwards to the
    /// base `intersection_size` implementation.
    ///
    /// @tparam           R The type of the first range.
    /// @tparam           S The type of the second range.
    /// @tparam ExecutionPolicy The type of the parallel execution policy.
    /// @tparam          _0 SFINAE to disambiguate from other 3 argument versions.
    ///
    /// @param           ep The parallel execution policy.
    /// @param            i The first range.
    /// @param            j The second range.
    ///
    /// @returns            The size of the intersected set.
    template <
      class R, class S, class ExecutionPolicy,
      std::enable_if_t<std::is_execution_policy_v<std::decay_t<ExecutionPolicy>>, void**> = nullptr>
    std::size_t intersection_size(ExecutionPolicy&& ep, R&& i, S&& j) {
      return intersection_size(std::forward<ExecutionPolicy>(ep), i.begin(), i.end(), j.begin(),
                               j.end());
    }

    // With custom comp function
    template <
      class ExecutionPolicy, class R, class S, class Compare,
      std::enable_if_t<std::is_execution_policy_v<std::decay_t<ExecutionPolicy>>, void**> = nullptr>
    std::size_t intersection_size(ExecutionPolicy&& ep, R&& i, S&& j, Compare&& cmp) {
      return intersection_size(std::forward<ExecutionPolicy>(ep), i.begin(), i.end(), j.begin(),
                               j.end(), std::forward<Compare>(cmp));
    }

    /// A convenience overload for `intersection_size`.
    ///
    /// This overload takes two iterators defining the first range, and a second
    /// range, and forwards to the base `intersection_size` implementation.
    ///
    /// @tparam           It1 The type of the first iterator.
    /// @tparam       Range The type of the second range.
    /// @tparam ExecutionPolicy The type of the parallel execution policy.
    /// @tparam          _0 SFINAE to disambiguate from other 4 argument versions.
    ///
    /// @param            i The beginning of the first range.
    /// @param           ie The end of the first range.
    /// @param            j The second range.
    /// @param           ep The parallel execution policy.
    ///
    /// @returns            The size of the intersected set.
    template <
      std::forward_iterator It1, class Range, class ExecutionPolicy,
      std::enable_if_t<std::is_execution_policy_v<std::decay_t<ExecutionPolicy>>, void**> = nullptr>
    std::size_t intersection_size(It1&& i, It1&& ie, Range&& j, ExecutionPolicy&& ep) {
      return intersection_size(std::forward<ExecutionPolicy>(ep), std::forward<It1>(i),
                               std::forward<It1>(ie), j.begin(), j.end());
    }

    /// A convenience overload for `intersection_size`.
    ///
    /// This overload takes two ranges as begin/end iterator pairs, and forwards to
    /// the base `intersection_size` with a sequential execution policy.
    ///
    /// @tparam           It1 The type of the first iterator.
    /// @tparam           It2 The type of the second iterator.
    ///
    /// @param            i The beginning of the first range.
    /// @param           ie The end of the first range.
    /// @param            j The beginning of the second range.
    /// @param           je The end of the second range.
    ///
    /// @returns            The size of the intersected set.
    template <class It1, class It2>
    std::size_t intersection_size(It1&& i, It1&& ie, It2&& j, It2&& je) {
      return intersection_size(std::execution::seq, std::forward<It1>(i), std::forward<It1>(ie),
                               std::forward<It2>(j), std::forward<It2>(je));
    }

    // With custom comp function
    template <std::forward_iterator It1, std::forward_iterator It2, class Compare>
    std::size_t intersection_size(It1&& i, It1&& ie, It2&& j, It2&& je, Compare&& cmp) {
      return intersection_size(std::execution::seq, std::forward<It1>(i), std::forward<It1>(ie),
                               std::forward<It2>(j), std::forward<It2>(je),
                               std::forward<Compare>(cmp));
    }

    /// A convenience overload for `intersection_size`.
    ///
    /// This overload takes two ranges and forwards to the base
    /// `intersection_size` implementation with a sequential execution policy.
    ///
    /// @tparam           R The type of the first range.
    /// @tparam           S The type of the second range.
    ///
    /// @param            i The first range.
    /// @param            j The second range.
    ///
    /// @returns            The size of the intersected set.
    template <class R, class S>
    std::size_t intersection_size(R&& i, S&& j) {
      return intersection_size(std::execution::seq, i.begin(), i.end(), j.begin(), j.end());
    }

    // With custom comp function
    template <class R, class S, class Compare>
    std::size_t intersection_size(R&& i, S&& j, Compare&& cmp) {
      return intersection_size(std::execution::seq, i.begin(), i.end(), j.begin(), j.end(),
                               std::forward<Compare>(cmp));
    }

    /// A convenience overload for `intersection_size`.
    ///
    /// This overload takes two iterators defining the first range, and a second
    /// range, and forwards to the base `intersection_size` implementation with a
    /// sequential execution policy.
    ///
    /// @tparam           A The type of the first iterator.
    /// @tparam           B The type of the second iterator.
    /// @tparam       Range The type of the second range.
    /// @tparam          _0 SFINAE to disambiguate from other 3 argument versions.
    ///
    /// @param            i The beginning of the first range.
    /// @param           ie The end of the first range.
    /// @param            j The second range.
    ///
    /// @returns            The size of the intersected set.
    template <std::forward_iterator It1, class Range>
    std::size_t intersection_size(It1&& i, It1&& ie, Range&& j) {
      return intersection_size(std::execution::seq, std::forward<It1>(i), std::forward<It1>(ie),
                               j.begin(), j.end());
    }
  } // namespace graph
} // namespace nw

#endif // NW_GRAPH_INTERSECTION_SIZE_HPP
