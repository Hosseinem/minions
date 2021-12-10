// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides modmer.
 */

#pragma once

#include <seqan3/std/algorithm>
#include <deque>

#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>

namespace seqan3::detail
{
// ---------------------------------------------------------------------------------------------------------------------
// modmer_view class
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief The type returned by modmer.
 * \tparam urng1_t The type of the underlying range, must model std::ranges::forward_range, the reference type must
 *                 model std::totally_ordered. The typical use case is that the reference type is the result of
 *                 seqan3::kmer_hash.
 * \implements std::ranges::view
 * \ingroup search_views
 *
 *
 * \note Most members of this class are generated by std::ranges::view_interface which is not yet documented here.

 */
template <std::ranges::view urng1_t, bool measure_distance = false>
class modmer_view : public std::ranges::view_interface<modmer_view<urng1_t>>
{
private:
    static_assert(std::ranges::forward_range<urng1_t>, "The modmer_view only works on forward_ranges.");
    static_assert(std::totally_ordered<std::ranges::range_reference_t<urng1_t>>,
                  "The reference type of the underlying range must model std::totally_ordered.");

    //!\brief Whether the given ranges are const_iterable
    static constexpr bool const_iterable = seqan3::const_iterable_range<urng1_t>;

    //!\brief The first underlying range.
    urng1_t urange1{};

    //!\brief The number of values in one window.
    size_t mod_used{};

    template <bool const_range>
    class basic_iterator;

    //!\brief The sentinel type of the modmer_view.
    using sentinel = std::default_sentinel_t;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    modmer_view()
        requires std::default_initializable<urng1_t>
        = default; //!< Defaulted.
    modmer_view(modmer_view const & rhs) = default; //!< Defaulted.
    modmer_view(modmer_view && rhs) = default; //!< Defaulted.
    modmer_view & operator=(modmer_view const & rhs) = default; //!< Defaulted.
    modmer_view & operator=(modmer_view && rhs) = default; //!< Defaulted.
    ~modmer_view() = default; //!< Defaulted.

    /*!\brief Construct from a view and a given number of values in one window.
    * \param[in] urange1     The input range to process. Must model std::ranges::viewable_range and
    *                        std::ranges::forward_range.
    * \param[in] mod_used The number of values in one window.
    */
    modmer_view(urng1_t urange1, size_t const mod_used) :
        urange1{std::move(urange1)},
        mod_used{mod_used}
    {}

    /*!\brief Construct from a non-view that can be view-wrapped and a given number of values in one window.
    * \tparam other_urng1_t  The type of another urange. Must model std::ranges::viewable_range and be constructible
                             from urng1_t.
    * \param[in] urange1     The input range to process. Must model std::ranges::viewable_range and
    *                        std::ranges::forward_range.
    * \param[in] mod_used The number of values in one window.
    */
    template <typename other_urng1_t>
    //!\cond
        requires (std::ranges::viewable_range<other_urng1_t> &&
                  std::constructible_from<urng1_t, ranges::ref_view<std::remove_reference_t<other_urng1_t>>>)
    //!\endcond
    modmer_view(other_urng1_t && urange1, size_t const mod_used) :
        urange1{std::views::all(std::forward<other_urng1_t>(urange1))},
        mod_used{mod_used}
    {}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the range.
     * \returns Iterator to the first element.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Strong exception guarantee.
     */
    basic_iterator<false> begin()
    {
        return {std::ranges::begin(urange1),
                std::ranges::end(urange1),
                mod_used};
    }

    //!\copydoc begin()
    basic_iterator<true> begin() const
    //!\cond
        requires const_iterable
    //!\endcond
    {
        return {std::ranges::cbegin(urange1),
                std::ranges::cend(urange1),
                mod_used};
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
     *
     * \details
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    sentinel end() const
    {
        return {};
    }
    //!\}
};

//!\brief Iterator for calculating modmers.
template <std::ranges::view urng1_t, bool measure_distance>
template <bool const_range>
class modmer_view<urng1_t, measure_distance>::basic_iterator
{
private:
    //!\brief The sentinel type of the first underlying range.
    using urng1_sentinel_t = maybe_const_sentinel_t<const_range, urng1_t>;
    //!\brief The iterator type of the first underlying range.
    using urng1_iterator_t = maybe_const_iterator_t<const_range, urng1_t>;

    template <bool>
    friend class basic_iterator;

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Type for distances between iterators.
    using difference_type = std::ranges::range_difference_t<urng1_t>;
    //!\brief Value type of this iterator.
    using value_type = std::ranges::range_value_t<urng1_t>;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Reference to `value_type`.
    using reference = value_type;
    //!\brief Tag this class as a forward iterator.
    using iterator_category = std::forward_iterator_tag;
    //!\brief Tag this class as a forward iterator.
    using iterator_concept = iterator_category;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    basic_iterator() = default; //!< Defaulted.
    basic_iterator(basic_iterator const &) = default; //!< Defaulted.
    basic_iterator(basic_iterator &&) = default; //!< Defaulted.
    basic_iterator & operator=(basic_iterator const &) = default; //!< Defaulted.
    basic_iterator & operator=(basic_iterator &&) = default; //!< Defaulted.
    ~basic_iterator() = default; //!< Defaulted.

    //!\brief Allow iterator on a const range to be constructible from an iterator over a non-const range.
    basic_iterator(basic_iterator<!const_range> const & it)
    //!\cond
        requires const_range
    //!\endcond
        : modmer_value{std::move(it.modmer_value)},
          urng1_iterator{std::move(it.urng1_iterator)},
          urng1_sentinel{std::move(it.urng1_sentinel)}
    {}

    /*!\brief Construct from begin and end iterators of a given range over std::totally_ordered values, and the number
              of values per window.
    * \param[in] urng1_iterator Iterator pointing to the first position of the first std::totally_ordered range.
    * \param[in] urng1_sentinel Iterator pointing to the last position of the first std::totally_ordered range.
    * \param[in] mod_used The number of values in one window.
    *
    * \details
    *
    * Looks at the number of values per window in two ranges, returns the smallest between both as modmer and
    * shifts then by one to repeat this action. If a modmer in consecutive windows is the same, it is returned only
    * once.
    */
    basic_iterator(urng1_iterator_t urng1_iterator,
                   urng1_sentinel_t urng1_sentinel,
                   size_t mod_used) :
        urng1_iterator{std::move(urng1_iterator)},
        urng1_sentinel{std::move(urng1_sentinel)},
        mod{mod_used}
    {
        size_t size = std::ranges::distance(urng1_iterator, urng1_sentinel);
        mod_used = std::min<size_t>(mod_used, size);

        first_modmer();
    }
    //!\}

    //!\anchor basic_iterator_comparison
    //!\name Comparison operators
    //!\{

    //!\brief Compare to another basic_iterator.
    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return (lhs.urng1_iterator == rhs.urng1_iterator);
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator!=(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to the sentinel of the modmer_view.
    friend bool operator==(basic_iterator const & lhs, sentinel const &)
    {
        return lhs.urng1_iterator == lhs.urng1_sentinel;
    }

    //!\brief Compare to the sentinel of the modmer_view.
    friend bool operator==(sentinel const & lhs, basic_iterator const & rhs)
    {
        return rhs == lhs;
    }

    //!\brief Compare to the sentinel of the modmer_view.
    friend bool operator!=(sentinel const & lhs, basic_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to the sentinel of the modmer_view.
    friend bool operator!=(basic_iterator const & lhs, sentinel const & rhs)
    {
        return !(lhs == rhs);
    }
    //!\}

    //!\brief Pre-increment.
    basic_iterator & operator++() noexcept
    {
        next_unique_modmer();
        return *this;
    }

    //!\brief Post-increment.
    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        next_unique_modmer();
        return tmp;
    }

    //!\brief Return the modmer.
    value_type operator*() const noexcept
    {
        return modmer_value;
    }

private:
    //!\brief The modmer value.
    value_type modmer_value{};

    //!\brief Iterator to the rightmost value of one window.
    urng1_iterator_t urng1_iterator{};
    //!brief Iterator to last element in range.
    urng1_sentinel_t urng1_sentinel{};

    size_t mod{};

    size_t distance{1};

    //!\brief Advances the window to the next position.
    void advance()
    {
        distance++;
        ++urng1_iterator;
    }

    void first_modmer()
    {
        if (!next_modmer())
            next_unique_modmer();
    }

    //!\brief Increments iterator by 1.
    void next_unique_modmer()
    {
        while (1)
        {
            advance();
            if (next_modmer())
                break;
        }
    }

    /*!\brief Calculates the next modmer value.
     * \returns True, if new modmer is found or end is reached. Otherwise returns false.
     */
    bool next_modmer()
    {
        if (urng1_iterator == urng1_sentinel)
            return true;

        if (*urng1_iterator % mod == 0)
        {
            if constexpr (measure_distance)
            {
                modmer_value = distance - 1;
                distance = 0;
            }
            else
            {
                modmer_value = *urng1_iterator;
            }
            return true;
        }

        return false;
    }
};

//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng1_t>
modmer_view(rng1_t &&, size_t const mod_used) -> modmer_view<std::views::all_t<rng1_t>>;

template <std::ranges::viewable_range rng1_t, bool m1>
modmer_view(rng1_t &&, size_t const mod_used) -> modmer_view<std::views::all_t<rng1_t>, m1>;


// ---------------------------------------------------------------------------------------------------------------------
// modmer_fn (adaptor definition)
// ---------------------------------------------------------------------------------------------------------------------

//![adaptor_def]
//!\brief modmer's range adaptor object type (non-closure).
//!\ingroup search_views
struct modmer_fn
{
    //!\brief Store the number of values in one window and return a range adaptor closure object.
    constexpr auto operator()(size_t const mod_used) const
    {
        return adaptor_from_functor{*this, mod_used};
    }

    /*!\brief Call the view's constructor with two arguments: the underlying view and an integer indicating how many
     *        values one window contains.
     * \tparam urng1_t        The type of the input range to process. Must model std::ranges::viewable_range.
     * \param[in] urange1     The input range to process. Must model std::ranges::viewable_range and
     *                        std::ranges::forward_range.
     * \param[in] mod_used The number of values in one window.
     * \returns  A range of converted values.
     */
    template <std::ranges::range urng1_t>
    constexpr auto operator()(urng1_t && urange1, size_t const mod_used) const
    {
        static_assert(std::ranges::viewable_range<urng1_t>,
                      "The range parameter to views::modmer cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng1_t>,
                      "The range parameter to views::modmer must model std::ranges::forward_range.");

        if (mod_used == 1) // Would just return urange1 without any changes
            throw std::invalid_argument{"The chosen mod_used is not valid. "
                                        "Please choose a value greater than 1 or use two ranges."};

        return modmer_view{urange1, mod_used};
    }
};
//![adaptor_def]

} // namespace seqan3::detail

/*!\brief Computes modmers for a range of comparable values. A modmer is ...
 * \tparam urng_t The type of the first range being processed. See below for requirements. [template
 *                 parameter is omitted in pipe notation]
 * \param[in] urange1 The range being processed. [parameter is omitted in pipe notation]
 * \param[in] mod_used The number of values in one window.
 * \returns A range of std::totally_ordered where each value is ... See below for the
 *          properties of the returned range.
 * \ingroup search_views
 *
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
 * |----------------------------------|:----------------------------------:|:--------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                      |
 * | std::ranges::forward_range       | *required*                         | *preserved*                      |
 * | std::ranges::bidirectional_range |                                    | *lost*                           |
 * | std::ranges::random_access_range |                                    | *lost*                           |
 * | std::ranges::contiguous_range    |                                    | *lost*                           |
 * |                                  |                                    |                                  |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
 * | std::ranges::view                |                                    | *guaranteed*                     |
 * | std::ranges::sized_range         |                                    | *lost*                           |
 * | std::ranges::common_range        |                                    | *lost*                           |
 * | std::ranges::output_range        |                                    | *lost*                           |
 * | seqan3::const_iterable_range     |                                    | *preserved*                      |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   | std::totally_ordered               | std::totally_ordered             |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 */
inline constexpr auto modmer = seqan3::detail::modmer_fn{};
