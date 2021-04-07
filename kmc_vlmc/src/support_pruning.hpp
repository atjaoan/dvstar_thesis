#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <kmc_file.h>

#include "kmer.hpp"

template <int kmer_size>
void output_root(std::vector<size_t> &root_counts,
                 kmer_sorter<kmer_size> &sorter) {
  size_t count = std::accumulate(root_counts.begin(), root_counts.end(), 0.0);
  VLMCKmer root{0, count,
                std::array<size_t, 4>{root_counts[0], root_counts[1],
                                      root_counts[2], root_counts[3]}};

  root.output(std::cout);
  sorter.push(root);
}

template <int kmer_size>
void output_node(VLMCKmer &prev_kmer, int diff_pos,
                 std::vector<std::vector<size_t>> &counters,
                 kmer_sorter<kmer_size> &sorter) {
  auto &next_character_counts = counters[diff_pos + 1];

  auto prefix_kmer = VLMCKmer::create_prefix_kmer(
      prev_kmer, diff_pos + 1, counters[diff_pos][0],
      std::array<size_t, 4>{next_character_counts[1], next_character_counts[2],
                            next_character_counts[3],
                            next_character_counts[4]});

  //  prefix_kmer.output(std::cout);
  sorter.push(prefix_kmer);
}

bool include_kmer(int length, size_t count) {
  int min_count = 1;
  int max_depth = 14;

  return length < max_depth && count >= min_count;
}

template <int kmer_size>
void process_kmer(VLMCKmer &current_kmer, VLMCKmer &prev_kmer,
                  std::vector<std::vector<size_t>> &counters,
                  std::vector<size_t> &root_counts,
                  kmer_sorter<kmer_size> &sorter) {
  // This should check what differs, and output the kmers that
  // don't match this level, with appropriate counters.
  //  std::cout << "current kmer:\t" << current_kmer.to_string() << std::endl;
  //  std::cout << "prev kmer:\t" << prev_kmer.to_string() << std::endl;

  auto diff_pos =
      VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);

  //  std::cout << "diff pos: " << diff_pos << std::endl;

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
  VLMCTranslator kmer_api(kmer_size);
  VLMCKmer kmer(kmer_size, 0, {});
  VLMCKmer prev_kmer(kmer_size, 0, {});

  int alphabet_size = 4;
  std::vector<std::vector<size_t>> counters(
      kmer_size + 2, std::vector<size_t>(alphabet_size + 1));
  std::vector<size_t> root_counts{};

  uint64 counter;

  while (kmer_database.ReadNextKmer(kmer_api, counter)) {
    kmer = kmer_api.construct_vlmc_kmer();
    kmer.count = counter;
    process_kmer(kmer, prev_kmer, counters, root_counts, sorter);

    prev_kmer = kmer;
  }
  VLMCTranslator api_epsilon(kmer_size);
  api_epsilon.from_string("N");
  VLMCKmer epsilon = api_epsilon.construct_vlmc_kmer();
  process_kmer(epsilon, prev_kmer, counters, root_counts, sorter);

  output_root(root_counts, sorter);
}
