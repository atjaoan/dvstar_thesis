#include <gtest/gtest.h>

#include "../src/support_pruning.hpp"

#include <kmc_file.h>

#include "read_helper.hpp"

using namespace vlmc;

class SupportPruningTests : public ::testing::Test {
protected:
  void SetUp() override {}

  int alphabet_size = 4;

  KMersPerLevel<uint64> counters{alphabet_size + 1, 7};

  std::shared_ptr<KmerContainer<ReverseKMerComparator<31>>> sorter =
      std::make_shared<OutOfCoreKmerContainer<ReverseKMerComparator<31>>>();

  vector_type local_kmers =
      std::make_shared<InCoreKmerContainer<ReverseKMerComparator<31>>>();

  const std::function<bool(int, uint64)> include_node =
      [](int length, uint64 count) -> bool { return true; };

  const std::function<void(const VLMCKmer &, int, uint64,
                           std::vector<uint64> &)>
      output_func = [&](const VLMCKmer &prev_kmer, int diff_pos, uint64 sum,
                        std::vector<uint64> &next_counts) {
        output_node(prev_kmer, diff_pos, sum, next_counts, local_kmers);
      };
};

TEST_F(SupportPruningTests, ProcessKMerOneDiff) {
  auto kmer = create_kmer("AAAAC");
  kmer.count = 2;
  counters[5][1] = 2;

  auto next_kmer = create_kmer("AAAAT");
  next_kmer.count = 5;

  process_kmer(next_kmer, kmer, counters, output_func, include_node, 1);

  local_kmers->for_each(
      [&](auto &kmer_) { EXPECT_EQ(kmer_.count, kmer.count); });

  increase_counts(counters, next_kmer);

  EXPECT_EQ(counters[5][3], next_kmer.count);
}

TEST_F(SupportPruningTests, ProcessKMerBigDiff) {
  auto kmer = create_kmer("ATTTT");
  counters[5][3] = 20;
  counters[5][2] = 20;
  counters[5][1] = 20;
  counters[5][0] = 20;

  auto next_kmer = create_kmer("CAAAA");
  next_kmer.count = 5;

  process_kmer(next_kmer, kmer, counters, output_func, include_node, 1);
  local_kmers->for_each([&](auto &kmer_) { EXPECT_EQ(kmer_.count, 20 * 4); });
  increase_counts(counters, next_kmer);

  EXPECT_EQ(counters[5][0], next_kmer.count);
}

TEST_F(SupportPruningTests, SortAll4Mers) {
  std::vector<VLMCKmer> kmers = {
      create_kmer("ACGT"), create_kmer("ACG"), create_kmer("ACTA"),
      create_kmer("ACT"),  create_kmer("AC"),  create_kmer("A"),
  };

  for (auto &kmer : kmers) {
    sorter->push(kmer);
  }
  sorter->sort();

  std::vector<VLMCKmer> stxxl_sorted_kmers{};
  sorter->for_each([&](VLMCKmer &kmer) { stxxl_sorted_kmers.push_back(kmer); });

  std::sort(kmers.begin(), kmers.end(), ReverseKMerComparator<31>());

  std::cout << "output" << std::endl;
  for (int i = 0; i < kmers.size(); i++) {
    kmers[i].output(std::cout);
    EXPECT_EQ(stxxl_sorted_kmers[i].to_string(), kmers[i].to_string());
  }
}

TEST_F(SupportPruningTests, PrefixSortTest) {
  std::vector<VLMCKmer> start_kmers = {
      create_kmer("ACGTA"),
      create_kmer("ACTAG"),
      create_kmer("TTTTA"),
  };

  for (int j = 1; j < start_kmers.size(); j++) {
    process_kmer(start_kmers[j], start_kmers[j - 1], counters, output_func,
                 include_node, 1);
  }

  local_kmers->for_each([&](auto &kmer) { sorter->push(kmer); });

  sorter->sort();

  std::vector<std::string> sorted_kmers = {"ACGT", "ACT",  "ACG",
                                           "AC",   "ACTA", "A"};

  int i = 0;
  sorter->for_each([&](VLMCKmer &kmer) {
    EXPECT_EQ(kmer.to_string(), sorted_kmers[i]);
    ++i;
  });
}
