#include <gtest/gtest.h>

#include "../src/support_pruning.hpp"

#include <kmc_file.h>

class SupportPruningTests : public ::testing::Test {
protected:
  void SetUp() override {}

  VLMCKmer current_kmer{15};
  VLMCKmer prev_kmer{15};
};

TEST_F(SupportPruningTests, DifferingPosition) {
  current_kmer.from_string("AAAAAAATCATCAGT");
  prev_kmer.from_string("AAAAAAAGTCGTCAG");

  auto diff_pos =
      VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, 7);
}

TEST_F(SupportPruningTests, NoDifferingPosition) {
  current_kmer.from_string("AAAAAAATCATCAGT");
  prev_kmer.from_string("AAAAAAATCATCAGT");

  auto diff_pos =
      VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, -1);
}

TEST_F(SupportPruningTests, AllDifferingPosition) {
  current_kmer.from_string("AAAAAAATCATCAGT");
  prev_kmer.from_string("TTCTCTATCTCTCTC");

  auto diff_pos =
      VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, 0);
}

TEST_F(SupportPruningTests, LongDifferingPosition) {
  VLMCKmer current_kmer{45};
  VLMCKmer prev_kmer{45};
  current_kmer.from_string("TACTAGCTACGATCATGCATGCATGCATGCAAAAAAATCATCAGT");
  prev_kmer.from_string("TACTAGCTACGATCATGCATGCATGCATGCAAAAAAAGTCGTCAG");

  auto diff_pos =
      VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, 37);
}

void check_prefixes(VLMCKmer &current_kmer, std::string &kmer_string) {
  for (int i = kmer_string.size(); i >= 0; i--) {
    auto prefix = VLMCKmer::create_prefix_kmer(current_kmer, (uint32)i, 0, {});
    std::string expected_string = kmer_string.substr(0, i);
    EXPECT_EQ(expected_string, prefix.to_string()) << i;
  }
}

TEST_F(SupportPruningTests, CreatePrefixKmer) {
  VLMCKmer current_kmer{16};
  std::string kmer_string{"TACTAGCTACGATCAT"};
  current_kmer.from_string(kmer_string);
  check_prefixes(current_kmer, kmer_string);
}

TEST_F(SupportPruningTests, CreatePrefixKmerOdd) {
  VLMCKmer current_kmer{15};
  std::string kmer_string{"TCGTACGACTGACAA"};
  current_kmer.from_string(kmer_string);

  check_prefixes(current_kmer, kmer_string);
}

TEST_F(SupportPruningTests, CreatePrefixKmerLong) {
  VLMCKmer current_kmer{45};
  std::string kmer_string{"TACTAGCTACGATCATTACTAGCTACGATCATTACTAGCTACGAT"};
  current_kmer.from_string(kmer_string);

  check_prefixes(current_kmer, kmer_string);
}

TEST_F(SupportPruningTests, SortTest) {
  VLMCKmer from_kmer{14};
  from_kmer.from_string("ACTGCTGATCGATC");

  VLMCKmer to_kmer{14};
  to_kmer.from_string("ACTGCTGATCGATA");

  EXPECT_LT(to_kmer, from_kmer);
}

TEST_F(SupportPruningTests, SortTestZero) {
  VLMCKmer from_kmer{14};
  from_kmer.from_string("ACTGCTGATCGATC");

  VLMCKmer to_kmer{0};

  EXPECT_LT(to_kmer, from_kmer);
}

TEST_F(SupportPruningTests, SortTestDiffLength) {
  VLMCKmer from_kmer{14};
  from_kmer.from_string("ACTGCTGATCGATC");

  VLMCKmer to_kmer{12};
  to_kmer.from_string("ACTGCTGATCGA");

  EXPECT_LT(to_kmer, from_kmer);
}
