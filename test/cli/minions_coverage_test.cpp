#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("minions coverage");
    std::string expected
    {
        "minions-coverage\n"
        "================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minimiser)
{
    cli_test_result result = execute_app("minions coverage --method minimiser -k 19 -w 19 ", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, gapped_minimiser)
{
    cli_test_result result = execute_app("minions coverage --method minimiser -k 19 -w 19 --shape 524223", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minstrobe)
{
    cli_test_result result = execute_app("minions coverage --method minstrobe --w_min 3 --w_max 5", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, syncmer)
{
    cli_test_result result = execute_app("minions coverage --method syncmer -K 6 -S 3", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, opensyncmer)
{
    cli_test_result result = execute_app("minions coverage --method opensyncmer -K 6 -S 3", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}


TEST_F(cli_test, modmer)
{
    cli_test_result result = execute_app("minions coverage --method modmer -k 19 -w 2 ", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, wrong_method)
{
    cli_test_result result = execute_app("minions coverage --method submer -k 19", data("example1.fasta"));
    std::string expected
    {
        "Error. Incorrect command line input for coverage. Validation failed "
        "for option --method: Value submer is not one of [kmer,minimiser,modmer,minstrobe,syncmer].\n"
    };

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.err, expected);
    EXPECT_EQ(result.out, std::string{});
}
