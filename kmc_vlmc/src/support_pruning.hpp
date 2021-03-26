#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <kmc_file.h>

#include "kmer.hpp"

template <int kmer_size>
void output_node(VLMCKmer &prev_kmer, int diff_pos,
                 std::vector<std::vector<size_t>> &counters,
                 kmer_sorter<kmer_size> &sorter) {
  for (int i = diff_pos; i < kmer_size; i++) {
    auto &next_character_counts = counters[diff_pos + 1];

    auto prefix_kmer = VLMCKmer::create_prefix_kmer(
        prev_kmer, diff_pos + 1, counters[diff_pos][0],
        std::vector<size_t>(next_character_counts.begin() + 1,
                            next_character_counts.end()));

    sorter.push(prefix_kmer);
  }
}

bool include_kmer(int length, size_t count) {
  int min_count = 10;
  int max_depth = 15;

  return length < max_depth && count >= min_count;
}

template <int kmer_size>
void process_kmer(VLMCKmer &current_kmer, VLMCKmer &prev_kmer,
                  std::vector<std::vector<size_t>> &counters,
                  std::vector<size_t> &root_counts,
                  kmer_sorter<kmer_size> &sorter) {
  // This should check what differs, and output the kmers that
  // don't match this level, with appropriate counters.

  auto diff_pos =
      VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);

//  std::cout << current_kmer.to_string() << " ";
//  std::cout << diff_pos;
//  std::cout << std::endl;

  if (diff_pos == -1) {
    return;
  }

  // Save the counts of the k-mers that are to be outputted:
  for (int i = diff_pos; i <= kmer_size; i++) {
    auto char_idx = prev_kmer.char_pos(i);
    int offset_char_idx = char_idx + 1;
    counters[i][offset_char_idx] = counters[i][0];
  }

  // Output kmers, reverse order to get counts right
  for (int i = kmer_size; i >= diff_pos; i--) {
    if (include_kmer(i, counters[i][0])) {
      output_node(prev_kmer, i, counters, sorter);
    }
    // Save root child counts
    if (i == 0) {
      root_counts.push_back(counters[i][0]);
    }

    // Reset counters for outputted kmer
    counters[i][0] = 0;
    std::fill(counters[i + 1].begin() + 1, counters[i + 1].end(), 0);
  }

  // The counts should be increased by this kmers counts for characters shorter
  // than the differing position.
  for (int i = 0; i <= kmer_size; i++) {
    counters[i][0] += current_kmer.count;
  }
}

template <int kmer_size>
void support_pruning(CKMCFile &kmer_database, kmer_sorter<kmer_size> &sorter) {
  VLMCKmer kmer(kmer_size);
  VLMCKmer prev_kmer(kmer_size);

  int alphabet_size = 4;
  std::vector<std::vector<size_t>> counters(
      kmer_size + 2, std::vector<size_t>(alphabet_size + 1));
  std::vector<size_t> root_counts{};

  uint64 counter;
  std::string str;

  while (kmer_database.ReadNextKmer(kmer, counter)) {
    kmer.count = counter;
    process_kmer(kmer, prev_kmer, counters, root_counts, sorter);

    kmer.to_string(str);
    prev_kmer = kmer;
  }
}
