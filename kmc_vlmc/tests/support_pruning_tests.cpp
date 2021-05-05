#include <gtest/gtest.h>

#include "../src/support_pruning.hpp"

#include <kmc_file.h>

#include "read_helper.hpp"

class SupportPruningTests : public ::testing::Test {
protected:
  void SetUp() override {}

  int alphabet_size = 4;

  KMersPerLevel<size_t> counters{alphabet_size + 1, 7};
  ;
  OutOfCoreKmerContainer<ReverseKMerComparator<31>> sorter{};

  const std::function<bool(int, size_t)> include_node =
      [](int length, size_t count) -> bool { return true; };
};

TEST_F(SupportPruningTests, ProcessKMerOneDiff) {
  auto kmer = create_kmer("AAAAC");
  kmer.count = 2;
  counters[4][0] = 2;

  auto next_kmer = create_kmer("AAAAT");
  next_kmer.count = 5;

  process_kmer(next_kmer, kmer, counters, sorter, include_node, 1);

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

  process_kmer(next_kmer, kmer, counters, sorter, include_node, 1);

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
      create_kmer("ACT"),  create_kmer("AC"),  create_kmer("A"),

  };

  for (auto &kmer : kmers) {
    sorter.push(kmer);
  }
  sorter.sort();

  std::vector<VLMCKmer> stxxl_sorted_kmers{};
  sorter.for_each([&](VLMCKmer &kmer) { stxxl_sorted_kmers.push_back(kmer); });

  std::sort(kmers.begin(), kmers.end(), ReverseKMerComparator<10>());

  std::cout << "output" << std::endl;
  for (int i = 0; i < kmers.size(); i++) {
    kmers[i].output(std::cout);
    ASSERT_EQ(stxxl_sorted_kmers[i].to_string(), kmers[i].to_string());
  }
}

TEST_F(SupportPruningTests, PrefixSortTest) {
  std::vector<VLMCKmer> start_kmers = {
      create_kmer("ACGTA"),
      create_kmer("ACTAG"),
      create_kmer("TTTTA"),
      create_kmer("TTTTT"),
  };

  std::cout << "input" << std::endl;
  for (int j = 1; j < start_kmers.size(); j++) {
    process_kmer(start_kmers[j], start_kmers[j - 1], counters, sorter,
                 include_node, 1);
  }

  sorter.sort();

  std::vector<std::string> sorted_kmers = {"ACGT", "ACT",  "ACG",
                                           "AC",   "ACTA", "A"};

  std::cout << "output" << std::endl;
  int i = 0;
  sorter.for_each([&](VLMCKmer &kmer) {
    EXPECT_EQ(kmer.to_string(), sorted_kmers[i]);
    ++i;
  });
}

TEST_F(SupportPruningTests, AllPrefixes) {
  auto prefixes = get_all_string_prefixes(2);

  EXPECT_EQ(prefixes[0], "AA");
  EXPECT_EQ(prefixes[1], "AC");
  EXPECT_EQ(prefixes[2], "AG");
  EXPECT_EQ(prefixes[3], "AT");
  EXPECT_EQ(prefixes[4], "CA");
  EXPECT_EQ(prefixes[5], "CC");
  EXPECT_EQ(prefixes[6], "CG");
  EXPECT_EQ(prefixes[7], "CT");
  EXPECT_EQ(prefixes[8], "GA");
  EXPECT_EQ(prefixes[9], "GC");
  EXPECT_EQ(prefixes[10], "GG");
  EXPECT_EQ(prefixes[11], "GT");
  EXPECT_EQ(prefixes[12], "TA");
  EXPECT_EQ(prefixes[13], "TC");
  EXPECT_EQ(prefixes[14], "TG");
  EXPECT_EQ(prefixes[15], "TT");
}