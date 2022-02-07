#include <forward_list>
#include <list>
#include <type_traits>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <seqan3/test/expect_range_eq.hpp>

#include <gtest/gtest.h>

#include "../../lib/seqan3/test/unit/range/iterator_test_template.hpp"

#include "syncmer.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

inline static constexpr auto smer_view = seqan3::views::kmer_hash(seqan3::ungapped{2});
inline static constexpr auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{5});

inline static constexpr auto syncmer_view = seqan3::views::syncmer(kmer_view,2,5);

using iterator_type = std::ranges::iterator_t< decltype(std::declval<seqan3::dna4_vector&>()
                                               | smer_view
                                               | syncmer_view)>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{2})) vec = text | smer_view;
    result_t expected_range{1,6,10,9,6,8,1,6,11,15,15,12,2};

    decltype(syncmer(seqan3::views::kmer_hash(text, seqan3::ungapped{2}), kmer_view, 2,5)) test_range =
    syncmer(vec, kmer_view, 2, 5);
};

using test_types = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_types, );

template <typename T>
class syncmer_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const,
                                                std::forward_list<seqan3::dna4>,
                                                std::forward_list<seqan3::dna4> const>;
TYPED_TEST_SUITE(syncmer_view_properties_test, underlying_range_types, );

class syncmer_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAA"_dna4};
    result_t result1{0, 0};
    result_t result1_open{0, 0};

    std::vector<seqan3::dna4> too_short_text{"AC"_dna4};

    // ACGG CGGC, GGCG, GCGA, CGAC, GACG, ACGT, CGTT, GTTT, TTTA, TTAG
    // CCGT GCCG  CGCC  TCGC  GTCG  CGTC  ACGT  AACG  AAAC  TAAA  CTAA
    // ACGG CGGC cgcc GCGA CGAC cgtc ACGT aacg aaac taaa ctaa
    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    result_t result3_ungapped{105,422,609,111,447,764,1010}; // ACGG, GGCG, GCGA, GACG, TTTA, TTAG
    result_t result3_open{105,422,111,447,764};
    result_t result3_ungapped_stop{105,422,609}; // ACGG, GGCG, GCGA, GACG
    result_t result3_open_stop{105,422};
    result_t result3_ungapped_start{111,447,764,1010};        // For start at second A, TTTA, TTAG
    result_t result3_open_start{111,447,764};
};

template <typename adaptor_t>
void compare_types(adaptor_t v)
{
    EXPECT_TRUE(std::ranges::input_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::view<decltype(v)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v), size_t>));
}

TYPED_TEST(syncmer_view_properties_test, concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG

    auto v = text | smer_view | syncmer_view;
    compare_types(v);
}

TYPED_TEST(syncmer_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped{109,438,865,111,447,764,1010};      // GTCG, TCGA, GACG, TTTA, TTAG
    EXPECT_RANGE_EQ(ungapped, text | smer_view | syncmer_view);
}

TEST_F(syncmer_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | smer_view | syncmer_view);
    auto empty_view = too_short_text | smer_view | syncmer_view;
    EXPECT_TRUE(std::ranges::empty(empty_view));
    EXPECT_RANGE_EQ(result3_ungapped, text3 | smer_view | syncmer_view);

    auto v1 = text1 | smer_view;
    auto v1.2 = text1 | kmer_view;
    EXPECT_RANGE_EQ(result1_open, (seqan3::detail::syncmer_view<decltype(v1), decltype(v1.2), true>(v1, v1.2, 2, 5)));
    auto v2 = text3 | smer_view;
    auto v2.2 = text3 | kmer_view;
    EXPECT_RANGE_EQ(result3_open, (seqan3::detail::syncmer_view<decltype(v2), decltype(v2.2), true>(v2, v2.2, 2, 5)));
}

TEST_F(syncmer_test, combinability)
{
    auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_ungapped_stop, text3 | stop_at_t | smer_view | syncmer_view);

    auto v1 = text3 | stop_at_t | smer_view;
    auto v2 = text3 | stop_at_t | kmer_view;

    EXPECT_RANGE_EQ(result3_open_stop, (seqan3::detail::syncmer_view<decltype(v1), decltype(v2), true>(v1, v2, 2, 5)));

    auto start_at_a = std::views::drop(6);
    EXPECT_RANGE_EQ(result3_ungapped_start, (seqan3::detail::syncmer_view{text3 | start_at_a | smer_view, text3 | start_at_a | kmer_view, 2, 5}));

    auto v3 = text3 | start_at_a | smer_view;
    auto v3.2 = text3 | start_at_a | kmer_view;

    EXPECT_RANGE_EQ(result3_open_start, (seqan3::detail::syncmer_view<decltype(v3), decltype(v3.2), true>(v3, 2,5)));
}
