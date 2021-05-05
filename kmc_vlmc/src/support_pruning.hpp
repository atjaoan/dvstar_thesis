#pragma once

#include <bitset>
#include <filesystem>
#include <fstream>
#include <functional>
#include <future>

#include <iostream>
#include <map>
#include <mutex>
#include <string>

#include <kmc_file.h>

#include "kmer.hpp"
#include "kmer_container.hpp"
#include "kmers_per_level.hpp"

std::mutex sorter_mutex{};
std::mutex print_mutex{};

void output_root(std::vector<size_t> &root_counts, KmerContainer<> &sorter) {
  size_t count = std::accumulate(root_counts.begin(), root_counts.end(), 0.0);
  VLMCKmer root{0, count,
                std::array<size_t, 4>{root_counts[1], root_counts[2],
                                      root_counts[3], root_counts[4]}};

  //  root.output(std::cout);
  std::lock_guard lock{sorter_mutex};
  sorter.push(root);
}


VLMCKmer output_node(VLMCKmer &prev_kmer, int length,
                     KMersPerLevel<size_t> &counters,
                     KmerContainer<> &output_kmers) {
  auto &next_counts = counters[length];

  size_t count =
      std::accumulate(next_counts.begin() + 1, next_counts.end(), 0.0);
  auto prefix_kmer = VLMCKmer::create_prefix_kmer(
      prev_kmer, length, count,
      std::array<size_t, 4>{next_counts[1], next_counts[2], next_counts[3],
                            next_counts[4]});

  //  prefix_kmer.output(std::cout);
//    std::lock_guard lock{sorter_mutex};
  output_kmers.push(prefix_kmer);

  return prefix_kmer;
}

void increase_counts(KMersPerLevel<size_t> &counters, size_t count) {
  for (int i = 0; i < counters.size(); i++) {
    counters[i][0] += count;
  }
}


VLMCKmer process_kmer(VLMCKmer &current_kmer, VLMCKmer &prev_kmer,
                      KMersPerLevel<size_t> &counters, KmerContainer<> &sorter,
                      const std::function<bool(int, size_t)> &include_node,
                      int prefix_length) {
  // This should check what differs, and output the kmers that
  // don't match this level, with appropriate counters.
  VLMCKmer output_kmer{};

  int diff_pos = 0;
  if (current_kmer.length == 0) {
    // Output the longest prefix of all the kmers in the current set
    diff_pos = prefix_length - 1;
  } else {
    diff_pos = VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);

    if (diff_pos == -1) {
      return output_kmer;
    }
  }

  // Save the counts of the k-mers that are to be outputted:
  for (int i = diff_pos; i < prev_kmer.length; i++) {
    auto char_idx = prev_kmer.char_pos(i);
    int offset_char_idx = char_idx + 1;
    counters[i][offset_char_idx] = counters[i][0];
  }

  // Output kmers, reverse order to get counts right
  for (size_t i = counters.size() - 1; i >= diff_pos; i--) {
    // For pseudo-counts, add 4 to counter
    if (include_node(i + 1, counters[i][0] + 4)) {
      output_node(prev_kmer, i, counters, sorter);
    }
  }

  for (int i = prev_kmer.length - 1; i >= diff_pos; i--) {
    // Reset counters for outputted kmer
    counters[i][0] = 0;
    std::fill(counters[i + 1].begin() + 1, counters[i + 1].end(), 0);
  }

  // The counts should be increased by this kmer's count
  increase_counts(counters, current_kmer.count);

  return output_kmer;
}

template <int kmer_size>
VLMCKmer run_kmer_subset(std::atomic<int> &running_tasks, int actual_kmer_size,
                         std::vector<VLMCKmer> kmers, int prefix_length,
                         KmerContainer<> &sorter,
                         const std::function<bool(int, size_t)> &include_node,
                         bool end) {

  running_tasks++;
  InCoreKmerContainer<ReverseKMerComparator<31>> local_sorter;
  int alphabet_size = 4;
  // actual_kmer_size + 2 as 0 is used for the root, and end is used for
  // next-counts of the longest kmers.
  KMersPerLevel<size_t> counters{alphabet_size + 1, kmer_size + 1};

  increase_counts(counters, kmers[0].count);

  for (int i = 1; i < kmers.size(); i++) {
    process_kmer(kmers[i], kmers[i - 1], counters, local_sorter,
                            include_node, prefix_length);
  }
  VLMCKmer epsilon(0, 0, {});
  auto kmer =
      process_kmer(epsilon, kmers[kmers.size() - 1], counters,
                              local_sorter, include_node, prefix_length);

  local_sorter.for_each([&](VLMCKmer &kmer_) { sorter.push(kmer_); });

  if (end) {
    output_root(counters[0], sorter);
  }

  running_tasks--;
  return kmer;
}

template <int kmer_size>
VLMCKmer run_kmc_kmer_subset(
    std::atomic<int> &running_tasks, int actual_kmer_size,
    std::string &kmc_db_name, VLMCKmer prefix, KmerContainer<> &sorter,
    const std::function<bool(int, size_t)> &include_node, bool end) {
  running_tasks++;

  InCoreKmerContainer<ReverseKMerComparator<31>> local_sorter;

  CKMCFile kmer_database;
  auto status = kmer_database.OpenForListing(kmc_db_name);

  auto n_kmers = kmer_database.KmerCount();

  VLMCTranslator kmer_api(actual_kmer_size);
  VLMCKmer kmer(actual_kmer_size, 0, {});
  VLMCKmer prev_kmer(actual_kmer_size, 0, {});
  uint64 counter;

  // Iterate kmer database until we find first kmer with the given prefix
  while (kmer_database.ReadNextKmer(kmer_api, counter)) {
    kmer = kmer_api.construct_vlmc_kmer();

    int diff = VLMCKmer::get_first_differing_position(prefix, kmer);
    if (diff == -1) {
      // Found the first kmer with this prefix.
      break;
    }
  }
  prev_kmer = kmer;
  prev_kmer.count = counter;

  int alphabet_size = 4;
  // actual_kmer_size + 2 as 0 is used for the root, and end is used for
  // next-counts of the longest kmers.
  KMersPerLevel<size_t> counters{alphabet_size + 1, kmer_size + 1};

  increase_counts(counters, counter);

  while (kmer_database.ReadNextKmer(kmer_api, counter)) {
    kmer = kmer_api.construct_vlmc_kmer();
    kmer.count = counter;

    int diff = VLMCKmer::get_first_differing_position(prefix, kmer);
    if (diff != -1) {
      // Found last k-mer with the given prefix.
      break;
    }

    process_kmer(kmer, prev_kmer, counters, local_sorter,
                            include_node, prefix.length);
    prev_kmer = kmer;
  }

  VLMCKmer epsilon(0, 0, {});
  auto return_kmer = process_kmer(
      epsilon, prev_kmer, counters, local_sorter, include_node, prefix.length);

  kmer_database.Close();

  std::lock_guard lock{sorter_mutex};
  local_sorter.for_each([&](VLMCKmer &kmer_) { sorter.push(kmer_); });

  if (end) {
    output_root(counters[0], sorter);
  }

  running_tasks--;

  return return_kmer;
}

VLMCKmer create_kmer(const std::string &kmer_string) {
  VLMCTranslator kmer{static_cast<int>(kmer_string.size())};
  if (!kmer_string.empty()) {
    kmer.from_string(kmer_string);
  }

  return kmer.construct_vlmc_kmer();
}

std::vector<std::string> get_all_string_prefixes(size_t length) {
  // All permutations of A, C, G, T

  std::vector<char> alphabet{'A', 'C', 'G', 'T'};
  size_t alphabet_size = 4;
  size_t n_prefixes = std::pow(alphabet_size, length);

  std::vector<std::string> prefixes(n_prefixes, std::string(length, ' '));

  for (size_t n = 0; n < length; n++) {
    size_t n_char = 0;
    for (auto char_ : alphabet) {
      size_t step_size = n_prefixes / std::pow(alphabet_size, n + 1);

      size_t start = step_size * n_char;
      size_t stop = step_size * (n_char + 1);

      while (stop <= n_prefixes) {
        for (size_t i = start; i < stop; i++) {
          prefixes[i][n] = char_;
        }
        start = start + step_size * alphabet_size;
        stop = stop + step_size * alphabet_size;
      }

      n_char++;
    }
  }

  return prefixes;
}

std::vector<VLMCKmer> get_prefixes(int prefix_length) {
  std::vector<VLMCKmer> kmers{};

  for (auto &prefix : get_all_string_prefixes(prefix_length)) {
    kmers.push_back(create_kmer(prefix));
  }

  return kmers;
}

template <int kmer_size>
void support_pruning(std::string &kmc_db_name, KmerContainer<> &sorter,
                     int actual_kmer_size,
                     const std::function<bool(int, size_t)> &include_node) {
  VLMCTranslator kmer_api(actual_kmer_size);
  VLMCKmer kmer(actual_kmer_size, 0, {});
  VLMCKmer prev_kmer(actual_kmer_size, 0, {});

  int alphabet_size = 4;
  // actual_kmer_size + 2 as 0 is used for the root, and end is used for
  // next-counts of the longest kmers.

  uint64 counter;

  int prefix_length = kmer_size;
  std::vector<std::future<VLMCKmer>> prefix_runs{};

  std::atomic<int> running_tasks = 0;
  int n_threads = 12; // TODO Replace with c++ get all threads
                      //  kmer_database.ReadNextKmer(kmer_api, counter);

  auto prefixes = get_prefixes(2);
  for (int i = 0; i < prefixes.size(); i++) {
    while (running_tasks >= n_threads) {
    }
    prefix_runs.push_back(std::async(
        std::launch::async, run_kmc_kmer_subset<kmer_size>,
        std::ref(running_tasks), actual_kmer_size, std::ref(kmc_db_name),
        prefixes[i], std::ref(sorter), include_node, false));
  }

  std::vector<VLMCKmer> kmers{};
  std::cout << "Task kmers" << std::endl;
  for (auto &future : prefix_runs) {
    kmers.push_back(future.get());
    kmers[kmers.size() - 1].output(std::cout);
  }
  std::cout << "Task run" << std::endl;
  run_kmer_subset<kmer_size>(running_tasks, actual_kmer_size, kmers, 1,
                             std::ref(sorter), include_node, true);
}

// template <int kmer_size>
// void sequential_support_pruning(
//     std::string &kmc_db_name, kmer_sorter<kmer_size> &sorter,
//     int actual_kmer_size,
//     const std::function<bool(int, size_t)> &include_node) {
//   CKMCFile kmer_database;
//   auto status = kmer_database.OpenForListing(kmc_db_name);
//
//   VLMCTranslator kmer_api(actual_kmer_size);
//   VLMCKmer kmer(actual_kmer_size, 0, {});
//   VLMCKmer prev_kmer(actual_kmer_size, 0, {});
//
//   int alphabet_size = 4;
//   // actual_kmer_size + 2 as 0 is used for the root, and end is used for
//   // next-counts of the longest kmers.
//   std::vector<std::vector<size_t>> counters(
//       actual_kmer_size + 2, std::vector<size_t>(alphabet_size + 1));
//
//   uint64 counter;
//
//   kmer_database.ReadNextKmer(kmer_api, counter);
//
//   prev_kmer = kmer_api.construct_vlmc_kmer();
//   prev_kmer.count = counter;
//
//   increase_counts(counters, counter);
//
//   while (kmer_database.ReadNextKmer(kmer_api, counter)) {
//     kmer = kmer_api.construct_vlmc_kmer();
//     kmer.count = counter;
//     process_kmer(kmer, prev_kmer, counters, sorter, include_node, 1);
//
//     prev_kmer = kmer;
//   }
//
//   VLMCKmer epsilon(0, 0, {});
//   process_kmer(epsilon, prev_kmer, counters, sorter, include_node, 1);
//   output_root(counters[0], sorter);
//
//   kmer_database.Close();
// }
