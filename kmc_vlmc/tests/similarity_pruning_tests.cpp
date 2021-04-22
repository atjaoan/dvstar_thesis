#include <gtest/gtest.h>

#include "../src/kmer.hpp"
#include "../src/similarity_pruning.hpp"

#include <kmc_file.h>

#include "read_helper.hpp"

class SimilarityPruningTests : public ::testing::Test {
protected:
  void SetUp() override {
    kmers_per_level = std::vector<std::array<PstKmer, 4>>(7);
  }
  VLMCKmer create_kmer(std::string kmer_string) {
    VLMCTranslator kmer{static_cast<int>(kmer_string.size())};
    if (kmer_string.size() > 0) {
      kmer.from_string(kmer_string);
    }

    return kmer.construct_vlmc_kmer();
  }

  std::vector<std::array<PstKmer, 4>> kmers_per_level;

  std::function<bool(double)> keep_node = [](double delta) -> bool {
    return delta <= 3.9075;
  };

};

TEST_F(SimilarityPruningTests, KullbackLiebler) {
  auto parent = create_kmer("TTTTT");
  parent.count = 6;
  parent.next_symbol_counts = {3, 1, 1, 1};

  auto child = create_kmer("TTTTTA");
  child.count = 3;
  child.next_symbol_counts = {0, 1, 1, 1};

  std::vector<double> parent_probs{4.0 / 10, 2.0 / 10, 2.0 / 10, 2.0 / 10};
  std::vector<double> child_probs{1.0 / 7, 2.0 / 7, 2.0 / 7, 2.0 / 7};

  double expected_kl =
      (child_probs[0] * std::log(child_probs[0] / parent_probs[0])) +
      (child_probs[1] * std::log(child_probs[1] / parent_probs[1])) +
      (child_probs[2] * std::log(child_probs[2] / parent_probs[2])) +
      (child_probs[3] * std::log(child_probs[3] / parent_probs[3]));
  expected_kl *= 3.0;

  double real_kl = kl_divergence(child, parent);

  EXPECT_EQ(expected_kl, real_kl);
}

TEST_F(SimilarityPruningTests, SimilarityPruneSameLevel) {
  auto prev_kmer = create_kmer("TTTTTT");
  auto kmer = create_kmer("ATTTTT");

  similarity_prune(prev_kmer, kmer, kmers_per_level, std::cout, keep_node);

  ASSERT_TRUE(kmers_per_level[6][0].real_child);
  EXPECT_EQ(kmers_per_level[6][0].kmer.to_string(), "ATTTTT");
}

TEST_F(SimilarityPruningTests, SimilarityPruneParent) {
  auto a_kmer = create_kmer("ATTTTT");
  a_kmer.count = 7;
  a_kmer.next_symbol_counts = {4, 1, 1, 1};

  auto c_kmer = create_kmer("CTTTTT");
  c_kmer.count = 3;
  c_kmer.next_symbol_counts = {1, 2, 0, 0};

  auto g_kmer = create_kmer("GTTTTT");
  g_kmer.count = 4;
  g_kmer.next_symbol_counts = {2, 1, 0, 1};

  auto t_kmer = create_kmer("TTTTTT");
  t_kmer.count = 1;
  t_kmer.next_symbol_counts = {0, 1, 0, 0};
  auto kmer = create_kmer("TTTTT");

  kmer.count = 15;
  kmer.next_symbol_counts = {7, 3, 4, 1};

  kmers_per_level[6][0] = {a_kmer, false, true};
  kmers_per_level[6][1] = {c_kmer, false, true};
  kmers_per_level[6][2] = {g_kmer, false, true};
  kmers_per_level[6][3] = {t_kmer, true, true};

  std::stringstream out;

  bool has_children = process_parent(a_kmer, kmer, kmers_per_level, out, keep_node);

  EXPECT_TRUE(has_children);

  EXPECT_EQ(kmers_per_level[6][0].real_child, false);
  EXPECT_EQ(kmers_per_level[6][1].real_child, false);
  EXPECT_EQ(kmers_per_level[6][2].real_child, false);
  EXPECT_EQ(kmers_per_level[6][3].real_child, false);

  std::string out_string = out.str();
  std::istringstream stream{out_string};
  for (std::string line; std::getline(stream, line);) {
    auto [kmer, count, next_symbol_counts] = read_line(line);
    EXPECT_EQ(kmer, "TTTTTT");
    EXPECT_EQ(count, t_kmer.count);
    EXPECT_EQ(next_symbol_counts[0], 0);
    EXPECT_EQ(next_symbol_counts[1], 1);
    EXPECT_EQ(next_symbol_counts[2], 0);
    EXPECT_EQ(next_symbol_counts[3], 0);
  }
}
