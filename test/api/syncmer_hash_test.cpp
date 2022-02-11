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

#include "syncmer_hash.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

using iterator_type = std::ranges::iterator_t<decltype(std::declval<seqan3::dna4_vector&>()
                                                       | syncmer_hash(2, 5, seqan3::seed{0}))>;

static constexpr auto ungapped_view = syncmer_hash(2, 5, seqan3::seed{0});

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = false;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    result_t expected_range{105,422,609,111,447,764,1010};

    using test_range_t = decltype(text | ungapped_view);
    test_range_t test_range = text | ungapped_view;
};

using test_type = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

template <typename T>
class syncmer_hash_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;

TYPED_TEST_SUITE(syncmer_hash_view_properties_test, underlying_range_types, );

class syncmer_hash_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAA"_dna4};
    result_t result1{0, 0};
//    result_t result1_open{0, 0};

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4}; //  Kmers:    ACGGC CGGCG GGCGA GCGAC CGACG GACGT ACGTT CGTTT GTTTA TTTAG
                                                            //  Hashed:    105,  422,  664,  609,  390,  539,  111,  447,  764,  1010
    result_t result3_ungapped{105,422,609,111,447,764,1010}; // Syncmers: ACGGC CGGCG       GCGAC             ACGTT CGTTT GTTTA TTTAG
//    result_t result3_open{105,422,111,447,764};  //          Openyncmers: ACGGC CGGCG                         ACGTT CGTTT GTTTA
    result_t result3_ungapped_stop{105,422,609}; //         Syncmer stop: ACGGC CGGCG       GCGAC
//    result_t result3_open_stop{105,422}; //             Opensyncmer stop: ACGGC CGGCG
    result_t result3_ungapped_start{111,447,764,1010}; //  Syncmer start:                                     ACGTT CGTTT GTTTA TTTAG
//    result_t result3_open_start{111,447,764};  //      Opensyncmer start:                                     ACGTT CGTTT GTTTA
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

    auto v = seqan3::detail::syncmer_hash(text, 5, 2, seqan3::seed{0});
    compare_types(v);
}

TYPED_TEST(syncmer_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped{109,438,865,111,447,764,1010};
    EXPECT_RANGE_EQ(ungapped, seqan3::detail::syncmer_hash(text, 5, 2, seqan3::seed{0}));
}

TEST_F(syncmer_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | ungapped_view);
    EXPECT_RANGE_EQ(result3_ungapped, text3 | ungapped_view);

//    EXPECT_RANGE_EQ(result1_open, (seqan3::detail::syncmer_hash<decltype(text1), true>(text1, 2, 5)));

//    EXPECT_RANGE_EQ(result3_open, (seqan3::detail::syncmer_hash<decltype(text3), true>(text3, 2, 5)));
}

TEST_F(syncmer_test, combinability)
{
    auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_ungapped_stop, text3 | stop_at_t | ungapped_view);

    auto v1 = text3 | stop_at_t | smer_view;
    auto v2 = text3 | stop_at_t | kmer_view;

//    EXPECT_RANGE_EQ(result3_open_stop, (seqan3::detail::syncmer_view<decltype(text3 | stop_at_t), true>(text3 | stop_at_t, 2, 5)));

    auto start_at_a = std::views::drop(6);
    EXPECT_RANGE_EQ(result3_ungapped_start, text3 | start_at_a | ungapped_view);

    auto v3 = text3 | start_at_a | smer_view;
    auto v3_2 = text3 | start_at_a | kmer_view;

//    EXPECT_RANGE_EQ(result3_open_start, (seqan3::detail::syncmer_hash<decltype(text3 | start_at_a), true>(text3 | start_at_a, 2, 5)));
}
