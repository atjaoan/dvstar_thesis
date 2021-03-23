#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <string>

#include <kmc_file.h>

#include "kmer.hpp"

void output_node(VLMCKmer &prev_kmer, int diff_pos,
                 std::vector<std::vector<size_t>> &counters) {
  // Create new kmer that is diff_pos long, and shares all characters with
  // prev_kmer. Set count on new kmer to counters[diff_pos][0] Set child couts
  // on new kmer to counters[diff_pos + 1][1,2,3,4]

  // Write node out, to file or to stxxl external vector
  // Every kmer that is longer than the diff pos should be outputted (to what?):
  // * Output to page files based on the reverse of the kmer?
  // * -To vector (but don't want to keep things in memory)-
  // * How does the KMC sorting work? - looks like a lot of work converting
  // their functions to a more general usecase
  // * Write to file or to
  // https://stxxl.org/tags/master/classstxxl_1_1vector.html
  //   ^ will make it easy to use stxxl's external memory sort later.
}

bool include_kmer(int length, size_t count) {
  int min_count = 10;
  int max_depth = 15;

  return length < max_depth && count >= min_count;
}

void process_kmer(VLMCKmer &current_kmer, VLMCKmer &prev_kmer, int kmer_length,
                  std::vector<std::vector<size_t>> &counters,
                  std::vector<size_t> &root_counts) {
  // This should check what differs, and output the kmers that
  // don't match this level, with appropriate counters.

  auto diff_pos =
      VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);

  std::cout << current_kmer.to_string() << " ";
  std::cout << diff_pos;
  std::cout << std::endl;

  // Save the counts of the k-mers that are to be outputted:
  for (int i = diff_pos; i <= kmer_length; i++) {
    auto char_idx = prev_kmer.char_pos(i);
    int offset_char_idx = char_idx + 1;
    counters[i][offset_char_idx] = counters[i][0];
  }

  // Output kmers, reverse order to get counts right
  for (int i = kmer_length; i >= diff_pos; i--) {
    if (include_kmer(i, counters[i][0])) {
      output_node(prev_kmer, i, counters);
    }
    // Save root child counts
    if (i == 0) {
      root_counts.push_back(counters[i][0]);
    }

    // Reset counters for outputted kmer
    counters[i][0] = 0;
    std::fill(counters[i + 1].begin() + 1, counters[i + 1].end(), 0);
  }

  // The count should be reset for every longer kmer than the diff

  // The counts should be increased by this kmers counts for characters shorter
  // than the differing position.

  for (int i = 0; i <= kmer_length; i++) {
    counters[i][0] += current_kmer.count;
  }
}

void support_pruning(CKMCFile &kmer_database, int kmer_size) {
  VLMCKmer kmer(kmer_size);
  VLMCKmer prev_kmer(kmer_size);

  int alphabet_size = 4;
  std::vector<std::vector<size_t>> counters(
      kmer_size + 1, std::vector<size_t>(alphabet_size + 1));
  std::vector<size_t> root_counts{};

  uint64 counter;
  std::string str;

  while (kmer_database.ReadNextKmer(kmer, counter)) {
    kmer.count = counter;
    process_kmer(kmer, prev_kmer, kmer_size, counters, root_counts);

    kmer.to_string(str);
    prev_kmer = kmer;
  }
  // Output last k-mer.
  // Output root node.
}
