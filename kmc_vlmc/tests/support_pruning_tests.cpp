#include <gtest/gtest.h>

#include "../src/support_pruning.hpp"

#include <kmc_file.h>

class SupportPruningTests : public ::testing::Test {
protected:
  void SetUp() override {}

  CKmerAPI current_kmer{15};
  CKmerAPI prev_kmer{15};
};

TEST_F(SupportPruningTests, DifferingPosition) {
  current_kmer.from_string("AAAAAAATCATCAGT");
  prev_kmer.from_string("AAAAAAAGTCGTCAG");

  auto diff_pos = get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, 7);
}

TEST_F(SupportPruningTests, NoDifferingPosition) {
  current_kmer.from_string("AAAAAAATCATCAGT");
  prev_kmer.from_string("AAAAAAATCATCAGT");

  auto diff_pos = get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, -1);
}

TEST_F(SupportPruningTests, AllDifferingPosition) {
  current_kmer.from_string("AAAAAAATCATCAGT");
  prev_kmer.from_string("TTCTCTATCTCTCTC");

  auto diff_pos = get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, 0);
}

TEST_F(SupportPruningTests, LongDifferingPosition) {
  CKmerAPI current_kmer{45};
  CKmerAPI prev_kmer{45};
  current_kmer.from_string("TACTAGCTACGATCATGCATGCATGCATGCAAAAAAATCATCAGT");
  prev_kmer.from_string("TACTAGCTACGATCATGCATGCATGCATGCAAAAAAAGTCGTCAG");

  auto diff_pos = get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, 37);
}