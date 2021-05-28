#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

#include <kmc_file.h>

#include "read_helper.hpp"

using namespace vlmc;

// TODO: This test checks for all k-mers in the sequence.
// This is not the same the count of all characters, as this also includes
// the last k "end"-mers, which give plus 1 to a couple of values, and plus
// k to the root node.

class CountCorrectnessTests : public ::testing::Test {
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

TEST_F(CountCorrectnessTests, CorrectCounts) {
  std::filesystem::path fasta_path{"../python-prototype/NC_022098.1.fasta"};
  std::ifstream fasta_stream{fasta_path};

  std::string sequence{};
  for (std::string line; std::getline(fasta_stream, line);) {
    if (line.empty()) {
      continue;
    }
    if (line[0] == '>') {
      continue;
    }

    sequence += line;
  }

  auto counts = get_kmer_counts(sequence);

  std::filesystem::path run_path{"out.txt"};
  std::ifstream run_stream{run_path};

  for (std::string line; std::getline(run_stream, line);) {
    if (line.empty()) {
      continue;
    }

    auto [kmer, count, next_symbol_counts] = read_line(line);

    EXPECT_EQ(counts[kmer], count) << kmer;
    EXPECT_EQ(counts[kmer + 'A'], next_symbol_counts[0]) << kmer + 'A';
    EXPECT_EQ(counts[kmer + 'C'], next_symbol_counts[1]) << kmer + 'C';
    EXPECT_EQ(counts[kmer + 'G'], next_symbol_counts[2]) << kmer + 'G';
    ASSERT_EQ(counts[kmer + 'T'], next_symbol_counts[3]) << kmer + 'T';
  }
}