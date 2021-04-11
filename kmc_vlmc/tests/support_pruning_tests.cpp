#include <gtest/gtest.h>

#include "../src/support_pruning.hpp"

#include <kmc_file.h>

class SupportPruningTests : public ::testing::Test {
protected:
  void SetUp() override {
    counters = std::vector<std::vector<size_t>>(
        7, std::vector<size_t>(alphabet_size + 1));
  }
  VLMCKmer create_kmer(std::string kmer_string) {
    VLMCTranslator kmer{static_cast<int>(kmer_string.size())};
    if (kmer_string.size() > 0) {
      kmer.from_string(kmer_string);
    }
    VLMCKmer kmer_ = kmer.construct_vlmc_kmer();
    kmer_.count = 10;

    return kmer_;
  }

  int alphabet_size = 4;

  std::vector<std::vector<size_t>> counters;
  kmer_sorter<5> sorter{ReverseKMerComparator<5>(), 64 * 1024 * 1024};
};

TEST_F(SupportPruningTests, ProcessKMerOneDiff) {
  auto kmer = create_kmer("AAAAC");
  kmer.count = 2;
  counters[4][0] = 2;

  auto next_kmer = create_kmer("AAAAT");
  next_kmer.count = 5;

  process_kmer(next_kmer, kmer, counters, sorter);

  EXPECT_EQ(counters[4][2], kmer.count);

  EXPECT_EQ(counters[0][0], next_kmer.count);
  EXPECT_EQ(counters[1][0], next_kmer.count);
  EXPECT_EQ(counters[2][0], next_kmer.count);
  EXPECT_EQ(counters[3][0], next_kmer.count);
  EXPECT_EQ(counters[4][0], next_kmer.count);
}

TEST_F(SupportPruningTests, ProcessKMerBigDiff) {
  auto kmer = create_kmer("ATTTT");
  counters[0][0] = 20;
  counters[1][0] = 20;
  counters[2][0] = 20;

  auto next_kmer = create_kmer("CAAAA");
  next_kmer.count = 5;

  process_kmer(next_kmer, kmer, counters, sorter);

  EXPECT_EQ(counters[0][1], 20);
  EXPECT_EQ(counters[1][4], 0);
  EXPECT_EQ(counters[2][4], 0);
  EXPECT_EQ(counters[3][4], 0);

  EXPECT_EQ(counters[0][0], next_kmer.count);
  EXPECT_EQ(counters[1][0], next_kmer.count);
  EXPECT_EQ(counters[2][0], next_kmer.count);
  EXPECT_EQ(counters[3][0], next_kmer.count);
  EXPECT_EQ(counters[4][0], next_kmer.count);
}

TEST_F(SupportPruningTests, SortAll4Mers) {
  std::vector<VLMCKmer> kmers = {
      create_kmer("ACGT"), create_kmer("ACG"), create_kmer("ACTA"),
      create_kmer("ACT"), create_kmer("AC"),  create_kmer("A"),

  };

  for (auto &kmer : kmers) {
    sorter.push(kmer);
  }
  sorter.sort();

  std::vector<VLMCKmer> stxxl_sorted_kmers{};
  while (!sorter.empty()) {
    stxxl_sorted_kmers.push_back(*sorter);
    ++sorter;
  }

  std::sort(kmers.begin(), kmers.end(), ReverseKMerComparator<10>());

  std::cout << "output" << std::endl;
  for (int i = 0; i < kmers.size(); i++) {
    kmers[i].output(std::cout);
    ASSERT_EQ(stxxl_sorted_kmers[i].to_string(), kmers[i].to_string());
  }
}

TEST_F(SupportPruningTests, PrefixSortTest) {
  std::vector<VLMCKmer> start_kmers = {
      create_kmer("AAAAA"), create_kmer("ACGTA"), create_kmer("ACTAG"),
      create_kmer("TTTTA"), create_kmer("TTTTT"),
  };

  std::cout << "input" << std::endl;
  for (int j = 1; j < start_kmers.size(); j++) {
    process_kmer(start_kmers[j], start_kmers[j - 1], counters, sorter);
  }

  sorter.sort();

  std::vector<std::string> sorted_kmers = {
      "ACGT", "ACT", "ACG", "AC", "ACTA", "A"
  };

  std::cout << "output" << std::endl;
  int i = 0;
  while (!sorter.empty()) {
    VLMCKmer kmer = *sorter;
    EXPECT_EQ(kmer.to_string(), sorted_kmers[i]);

    ++i;
    ++sorter;
  }
}
