#include <gtest/gtest.h>

#include "../src/support_pruning.hpp"

#include <kmc_file.h>

class SupportPruningTests : public ::testing::Test {
protected:
  void SetUp() override {}

  VLMCKmer create_kmer(std::string kmer_string) {
    VLMCTranslator kmer{kmer_string.size()};
    if (kmer_string.size() > 0) {
      kmer.from_string(kmer_string);
    }

    return kmer.construct_vlmc_kmer();
  }

  VLMCKmer current_kmer{15, 0, {}};
  VLMCKmer prev_kmer{15, 0, {}};
};

TEST_F(SupportPruningTests, VLMCKmerConstructor) {
  std::string kmer_string{"ACTGCTGATCGA"};
  auto created = create_kmer(kmer_string);

  EXPECT_EQ(created.to_string(), kmer_string);
}

TEST_F(SupportPruningTests, VLMCKmerConstructorOdd) {
  std::string kmer_string{"ACTGCTGATCGAC"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(created.to_string(), kmer_string);
}

TEST_F(SupportPruningTests, VLMCKmerConstructor11) {
  std::string kmer_string{"ACTGCTGATCG"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(created.to_string(), kmer_string);
}

TEST_F(SupportPruningTests, VLMCKmerConstructor10) {
  std::string kmer_string{"ACTGCTGATC"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(created.to_string(), kmer_string);
}

TEST_F(SupportPruningTests, VLMCKmerConstructor9) {
  std::string kmer_string{"ACTGCTGAT"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(created.to_string(), kmer_string);
}

TEST_F(SupportPruningTests, IsPod) {
  EXPECT_FALSE(std::is_pod<CKmerAPI>::value);
  EXPECT_TRUE(std::is_pod<VLMCKmer>::value);
}

TEST_F(SupportPruningTests, DifferingPosition) {
  auto current_kmer = create_kmer("AAAAAAATCATCAGT");
  auto prev_kmer = create_kmer("AAAAAAAGTCGTCAG");

  auto diff_pos =
      VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, 7);
}

TEST_F(SupportPruningTests, NoDifferingPosition) {
  auto current_kmer = create_kmer("AAAAAAATCATCAGT");
  auto prev_kmer = create_kmer("AAAAAAATCATCAGT");

  auto diff_pos =
      VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, -1);
}

TEST_F(SupportPruningTests, AllDifferingPosition) {
  auto current_kmer = create_kmer("AAAAAAATCATCAGT");
  auto prev_kmer = create_kmer("TTCTCTATCTCTCTC");

  auto diff_pos =
      VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);
  EXPECT_EQ(diff_pos, 0);
}

TEST_F(SupportPruningTests, LongDifferingPosition) {
  auto current_kmer =
      create_kmer("TACTAGCTACGATCATGCATGCATGCATGCAAAAAAATCATCAGT");
  auto prev_kmer = create_kmer("TACTAGCTACGATCATGCATGCATGCATGCAAAAAAAGTCGTCAG");

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
  std::string kmer_string{"TACTAGCTACGATCAT"};
  auto current_kmer = create_kmer(kmer_string);

  check_prefixes(current_kmer, kmer_string);
}

TEST_F(SupportPruningTests, CreatePrefixKmerOdd) {
  std::string kmer_string{"TCGTACGACTGACAA"};
  auto current_kmer = create_kmer(kmer_string);

  check_prefixes(current_kmer, kmer_string);
}

TEST_F(SupportPruningTests, CreatePrefixKmerLong) {
  std::string kmer_string{"TACTAGCTACGATCATTACTAGCTACGATCATTACTAGCTACGAT"};
  auto current_kmer = create_kmer(kmer_string);

  check_prefixes(current_kmer, kmer_string);
}

TEST_F(SupportPruningTests, SortTest) {
  auto from_kmer = create_kmer("ACTGCTGATCGATC");
  auto to_kmer = create_kmer("ACTGCTGATCGATA");

  EXPECT_LT(to_kmer, from_kmer);
}

TEST_F(SupportPruningTests, SortTestZero) {
  auto from_kmer = create_kmer("ACTGCTGATCGATC");
  auto to_kmer = create_kmer("");

  EXPECT_LT(to_kmer, from_kmer);
}

TEST_F(SupportPruningTests, SortTestDiffLength) {
  auto from_kmer = create_kmer("ACTGCTGATCGATC");
  auto to_kmer = create_kmer("ACTGCTGATCGA");

  EXPECT_LT(to_kmer, from_kmer);
}
