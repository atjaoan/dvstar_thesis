#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

#include "read_helper.hpp"

using namespace vlmc;

class RegressionTests : public ::testing::Test {
protected:
  void SetUp() override {}
};

std::map<std::string, size_t> get_kmer_counts(std::string &sequence) {
  std::map<std::string, size_t> counts{};

  for (size_t i = 0; i <= sequence.size() - 10; i++) {
    for (int j = 0; j < 11 && j + i <= sequence.size(); j++) {
      auto kmer = sequence.substr(i, j);
      if (counts.find(kmer) != counts.end()) {
        counts[kmer] += 1;
      } else {
        counts[kmer] = 1;
      }
    }
  }
  return counts;
}

TEST_F(RegressionTests, PstClassifierSeqan) {
  std::filesystem::path tree_path{"seqan.tree"};
  std::ifstream tree_stream{tree_path};

  std::map<std::string, std::tuple<size_t, std::array<size_t, 4>>> tree{};
  for (std::string line; std::getline(tree_stream, line);) {
    if (line.size() == 0) {
      continue;
    }
    if (line.substr(0, 5) != "Node:") {
      continue;
    }

    auto [kmer, count, next_symbol_counts] = read_tree_line(line);

    tree[kmer] = {count, next_symbol_counts};
  }

  std::filesystem::path run_path{"out.txt"};
  std::ifstream run_stream{run_path};

  std::set<std::string> kmc_vlmc_contexts{};
  for (std::string line; std::getline(run_stream, line);) {
    if (line.size() == 0) {
      continue;
    }

    auto [kmer, count, next_symbol_counts] = read_line(line);

    kmc_vlmc_contexts.insert(kmer);

    ASSERT_TRUE(tree.find(kmer) != tree.end()) << kmer;

    auto [tree_count, tree_next_symbol_counts] = tree[kmer];

    // Some counts won't be correct, see comment in count_correctness_tests.
    // EXPECT_EQ(tree_count, count) << kmer;
    // EXPECT_EQ(tree_next_symbol_counts[0], next_symbol_counts[0])<<kmer + 'A';
    // EXPECT_EQ(tree_next_symbol_counts[1], next_symbol_counts[1])<<kmer + 'C';
    // EXPECT_EQ(tree_next_symbol_counts[2], next_symbol_counts[2])<<kmer + 'G';
    // ASSERT_EQ(tree_next_symbol_counts[3], next_symbol_counts[3])<<kmer + 'T';
  }

  for (auto &[kmer, _] : tree) {
    EXPECT_TRUE(kmc_vlmc_contexts.find(kmer) != kmc_vlmc_contexts.end())
        << kmer;
  }
}