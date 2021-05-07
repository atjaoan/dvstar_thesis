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

typedef std::vector<VLMCKmer> vector_type;

void output_root(std::vector<size_t> &root_counts,
                 std::shared_ptr<KmerContainer<>> &container) {
  size_t count = std::accumulate(root_counts.begin(), root_counts.end(), 0.0);
  VLMCKmer root{0, count,
                std::array<uint64, 4>{root_counts[1], root_counts[2],
                                      root_counts[3], root_counts[4]}};

  //  root.output(std::cout);
  //  std::lock_guard lock{sorter_mutex};
  container->push(root);
}

VLMCKmer output_node(const VLMCKmer &prev_kmer, int diff_pos,
                     KMersPerLevel<size_t> &counters,
                     vector_type &output_kmers) {
  auto &next_counts = counters[diff_pos + 1];

  auto prefix_kmer = VLMCKmer::create_prefix_kmer(
      prev_kmer, diff_pos + 1, counters[diff_pos][0],
      std::array<uint64, 4>{next_counts[1], next_counts[2], next_counts[3],
                            next_counts[4]});

  //  prefix_kmer.output(std::cout);
  //    std::lock_guard lock{sorter_mutex};
  output_kmers.push_back(prefix_kmer);

  return prefix_kmer;
}

void increase_counts(KMersPerLevel<size_t> &counters, size_t count) {
  for (int i = 0; i < counters.size(); i++) {
    counters[i][0] += count;
  }
}

VLMCKmer process_kmer(const VLMCKmer &current_kmer, const VLMCKmer &prev_kmer,
                      KMersPerLevel<size_t> &counters,
                      vector_type &output_kmers,
                      const std::function<bool(int, size_t)> &include_node,
                      int prefix_length) {
  // This should check what differs, and output the kmers that
  // don't match this level, with appropriate counters.
  VLMCKmer output_kmer{0, 0, {}};

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
  for (int i = prev_kmer.length - 1; i >= diff_pos; i--) {
    // For pseudo-counts, add 4 to counter
    // And the length has to be at least one shorter than the prev_kmer
    // to allow variable k iterations
    if (i < prev_kmer.length - 1 && include_node(i + 1, counters[i][0] + 4)) {
      output_kmer = output_node(prev_kmer, i, counters, output_kmers);
    }

    // Reset counters for potentially outputted kmer
    counters[i][0] = 0;
    std::fill(counters[i + 1].begin() + 1, counters[i + 1].end(), 0);
  }

  return output_kmer;
}

template <int kmer_size>
std::tuple<VLMCKmer, vector_type>
run_kmer_subset(std::atomic<int> &running_tasks, int actual_kmer_size,
                std::shared_ptr<KmerContainer<KMerComparator<31>>> &in_kmers,
                int prefix_length, std::shared_ptr<KmerContainer<>> &container,
                const std::function<bool(int, size_t)> &include_node, bool end,
                bool sorted = false) {
  if (in_kmers->size() == 0) {
    return {VLMCKmer{0, 0, {}}, {}};
  }
  if (!sorted) {
    in_kmers->sort();
  }

  vector_type local_kmers{};
  int alphabet_size = 4;
  KMersPerLevel<size_t> counters{alphabet_size + 1, kmer_size + 1};

  VLMCKmer prev_kmer{0, 0, {}};
  in_kmers->for_each([&](VLMCKmer &kmer) {
    if (prev_kmer.length != 0) {
      process_kmer(kmer, prev_kmer, counters, local_kmers, include_node,
                   prefix_length);
    }

    // The counts should be increased by this kmer's count
    increase_counts(counters, kmer.count);

    prev_kmer = kmer;
  });

  VLMCKmer epsilon(0, 0, {});
  auto kmer = process_kmer(epsilon, prev_kmer, counters, local_kmers,
                           include_node, prefix_length);

  in_kmers->clear();

  //  std::lock_guard lock{sorter_mutex};
  //  for (auto &kmer_ : local_kmers) {
  //    //    kmer_.output(std::cout);
  //    container->push(kmer_);
  //  }

  if (end) {
    output_root(counters[0], container);
  }

  running_tasks--;
  return {kmer, local_kmers};
}

template <int kmer_size>
void support_pruning(CKMCFile &kmer_database,
                     std::shared_ptr<KmerContainer<>> &container,
                     int actual_kmer_size,
                     const std::function<bool(int, size_t)> &include_node,
                     const std::string &in_or_out_of_core) {
  VLMCTranslator kmer_api(actual_kmer_size);
  uint64 counter;

  VLMCKmer kmer(actual_kmer_size, 0, {});

  int prefix_length = 2; // TODO determine this based on task size in some way
  int n_tasks = std::pow(4, prefix_length);

  std::vector<std::shared_ptr<KmerContainer<KMerComparator<31>>>> tasks(
      n_tasks);
  for (int i = 0; i < n_tasks; i++) {
    if (in_or_out_of_core == "internal") {
      tasks[i] = std::make_shared<InCoreKmerContainer<KMerComparator<31>>>(
          Iteration::sequential);
    } else {
      tasks[i] = std::make_shared<OutOfCoreKmerContainer<KMerComparator<31>>>();
    }
  }

  while (kmer_database.ReadNextKmer(kmer_api, counter)) {
    kmer = kmer_api.construct_vlmc_kmer();
    kmer.count = counter;

    size_t prefix_idx = kmer.get_prefix_index(prefix_length);

    tasks[prefix_idx]->push(kmer);

    //    std::cout << kmer.to_string() << " " << prefix_idx << std::endl;
  }

  std::atomic<int> running_tasks = 0;
  int n_threads = 12; // TODO Replace with c++ get all threads

  if (in_or_out_of_core == "external") {
    // In the exernal case, we don't have a sequential sort, so it's
    // probably better to sort before creating parallel tasks.
    for (auto &kmers : tasks) {
      while (running_tasks >= n_threads) {
      }
      if (kmers->size() == 0) {
        continue;
      }
      kmers->sort();
    }
  }

  std::vector<std::future<std::tuple<VLMCKmer, vector_type>>> prefix_runs{};
  for (auto &kmers : tasks) {
    while (running_tasks >= n_threads) {
    }
    if (kmers->size() == 0) {
      continue;
    }
    running_tasks++;
    prefix_runs.push_back(std::async(
        run_kmer_subset<kmer_size>, std::ref(running_tasks), actual_kmer_size,
        std::ref(kmers), prefix_length, std::ref(container), include_node,
        false, in_or_out_of_core == "external"));
  }

  std::shared_ptr<KmerContainer<KMerComparator<31>>> merge_kmers =
      std::make_shared<InCoreKmerContainer<KMerComparator<31>>>();

  for (auto &future : prefix_runs) {
    auto [kmer_, kmers] = future.get();

    for (auto &kmer_ : kmers) {
      //    kmer_.output(std::cout);
      container->push(kmer_);
    }
    if (kmer_.length != 0) {
      merge_kmers->push(kmer_);
    }
  }

  auto [kmer_, kmers] =
      run_kmer_subset<kmer_size>(running_tasks, actual_kmer_size, merge_kmers,
                                 1, std::ref(container), include_node, true);
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
