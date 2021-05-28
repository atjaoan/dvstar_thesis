#pragma once

#include <atomic>
#include <bitset>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <string>
#include <tuple>

#include <kmc_file.h>

#include "kmer.hpp"
#include "kmer_container.hpp"
#include "kmers_per_level.hpp"

namespace vlmc {

using vector_type = std::shared_ptr<KmerContainer<>>;

void output_root(std::vector<uint64> &root_counts, vector_type &container) {
  uint64 count = std::accumulate(root_counts.begin(), root_counts.end(), 0);
  VLMCKmer root{0, count,
                std::array<uint64, 4>{root_counts[1], root_counts[2],
                                      root_counts[3], root_counts[4]}};

  //  root.output(std::cout);
  //  std::lock_guard lock{sorter_mutex};
  container->push(root);
}

void output_node(const VLMCKmer &prev_kmer, int diff_pos, uint64 sum,
                 std::vector<uint64> &next_counts, vector_type &output_kmers) {
  auto prefix_kmer = VLMCKmer::create_prefix_kmer(
      prev_kmer, diff_pos + 1, sum,
      std::array<uint64, 4>{next_counts[1], next_counts[2], next_counts[3],
                            next_counts[4]});

  //  prefix_kmer.output(std::cout);
  //    std::lock_guard lock{sorter_mutex};
  output_kmers->push(prefix_kmer);
}

void increase_counts(KMersPerLevel<uint64> &counters, const VLMCKmer &kmer) {
  uchar char_idx = kmer.char_pos(kmer.length - 1);
  counters[kmer.length][char_idx + 1] += kmer.count;
}

void process_kmer(const VLMCKmer &current_kmer, const VLMCKmer &prev_kmer,
                  KMersPerLevel<uint64> &counters,
                  const std::function<void(const VLMCKmer &, int, uint64,
                                           std::vector<uint64> &)> &output_func,
                  const std::function<bool(int, size_t)> &include_node,
                  int prefix_length) {
  // This should check what differs, and output the kmers that
  // don't match this level, with appropriate counters.

  int32 diff_pos = 0;
  if (current_kmer.length == 0) {
    // Output the longest prefix of all the kmers in the current set
    diff_pos = prefix_length - 1;
  } else {
    diff_pos = VLMCKmer::get_first_differing_position(current_kmer, prev_kmer);

    if (diff_pos == -1) {
      return;
    }
  }

  // Output kmers, longest to shortest to get counts right
  for (int i = prev_kmer.length - 1; i >= diff_pos; i--) {
    // For pseudo-counts, add 4 to counter
    // And the length has to be at least one shorter than the prev_kmer
    // to allow variable k iterations
    auto &next_counts = counters[i + 1];

    uint64 sum =
        next_counts[1] + next_counts[2] + next_counts[3] + next_counts[4];

    if (i < prev_kmer.length - 1 && include_node(i + 1, sum + 4)) {
      output_func(prev_kmer, i, sum, next_counts);
      //      output_node(prev_kmer, i, sum, next_counts, output_kmers);
    }

    auto char_idx = prev_kmer.char_pos(i);
    counters[i][char_idx + 1] = sum;

    // Reset counters for potentially outputted kmer
    counters[i][0] = 0;
    next_counts[1] = next_counts[2] = next_counts[3] = next_counts[4] = 0;
  }
}

template <int kmer_size>
std::tuple<VLMCKmer, vector_type>
run_kmer_subset(std::atomic<int> &running_tasks, int actual_kmer_size,
                std::shared_ptr<KmerContainer<KMerComparator<31>>> &in_kmers,
                int prefix_length, vector_type &container,
                const std::function<bool(int, size_t)> &include_node,
                bool end) {
  vector_type local_kmers = std::make_shared<InCoreKmerContainer<>>();

  if (in_kmers->size() == 0) {
    running_tasks--;
    return {VLMCKmer{0, 0, {}}, local_kmers};
  }

  auto output_func = [&](const VLMCKmer &prev_kmer, int diff_pos, uint64 sum,
                         std::vector<uint64> &next_counts) {
    output_node(prev_kmer, diff_pos, sum, next_counts, local_kmers);
  };

  int alphabet_size = 4;
  KMersPerLevel<uint64> counters{alphabet_size + 1, kmer_size + 1};

  VLMCKmer prev_kmer{0, 0, {}};
  in_kmers->for_each([&](VLMCKmer &kmer) {
    if (prev_kmer.length != 0) {
      process_kmer(kmer, prev_kmer, counters, output_func, include_node,
                   prefix_length);
    }

    // The counts should be increased by this kmer's count
    increase_counts(counters, kmer);

    prev_kmer = kmer;
  });

  VLMCKmer epsilon(0, 0, {});
  process_kmer(epsilon, prev_kmer, counters, output_func, include_node,
               prefix_length);
  in_kmers->clear();

  if (end) {
    output_root(counters[0], container);
  }

  running_tasks--;
  return {VLMCKmer{0, 0, {}}, local_kmers};
}

template <int kmer_size>
void support_pruning(CKMCFile &kmer_database, vector_type &container,
                     int actual_kmer_size,
                     const std::function<bool(int, size_t)> &include_node,
                     const std::string &in_or_out_of_core) {
  VLMCTranslator kmer_api(actual_kmer_size);
  uint64 counter;

  VLMCKmer kmer(actual_kmer_size, 0, {});

  int prefix_length = 3; // TODO determine this based on task size in some way
  int n_tasks = std::pow(4, prefix_length);

  std::vector<std::shared_ptr<KmerContainer<KMerComparator<31>>>> tasks(
      n_tasks);
  for (int i = 0; i < n_tasks; i++) {
    if (in_or_out_of_core == "internal") {
      tasks[i] = std::make_shared<InCoreKmerContainer<KMerComparator<31>>>(
          Iteration::sequential);
    } else {
      tasks[i] = std::make_shared<
          OutOfCoreSequentialKmerContainer<KMerComparator<31>>>();
    }
  }

  std::atomic<int> running_tasks = 0;
  int n_threads = 12; // TODO Replace with c++ get all threads
  std::vector<std::future<std::tuple<VLMCKmer, vector_type>>> prefix_runs{};

  size_t current_prefix_idx = 0;
  while (kmer_database.ReadNextKmer(kmer_api, counter)) {
    kmer_api.construct_vlmc_kmer(kmer);
    kmer.count = counter;

    size_t prefix_idx = kmer.get_prefix_index(prefix_length);
    // std::cout << kmer.to_string() << " " << prefix_idx << std::endl;

    if (current_prefix_idx != prefix_idx) {
      running_tasks++;
      prefix_runs.push_back(
          std::async(std::launch::async, run_kmer_subset<kmer_size>,
                     std::ref(running_tasks), actual_kmer_size,
                     std::ref(tasks[current_prefix_idx]), prefix_length,
                     std::ref(container), include_node, false));
    }

    tasks[prefix_idx]->push(kmer);
    current_prefix_idx = prefix_idx;
  }
  prefix_runs.push_back(std::async(
      std::launch::async, run_kmer_subset<kmer_size>, std::ref(running_tasks),
      actual_kmer_size, std::ref(tasks[current_prefix_idx]), prefix_length,
      std::ref(container), include_node, false));

  std::shared_ptr<KmerContainer<KMerComparator<31>>> merge_kmers =
      std::make_shared<InCoreKmerContainer<KMerComparator<31>>>();

  for (auto &future : prefix_runs) {
    auto [kmer_, kmers] = future.get();

    kmers->for_each([&](auto &kmer_) { container->push(kmer_); });

    if (kmer_.length != 0) {
      merge_kmers->push(kmer_);
    }
  }

  auto [kmer_, kmers] =
      run_kmer_subset<kmer_size>(running_tasks, actual_kmer_size, merge_kmers,
                                 1, std::ref(container), include_node, true);

  kmers->for_each([&](auto &kmer_) { container->push(kmer_); });
}

template <int kmer_size>
void sequential_support_pruning(
    CKMCFile &kmer_database, vector_type &container, int actual_kmer_size,
    const std::function<bool(int, size_t)> &include_node) {

  VLMCTranslator kmer_api(actual_kmer_size);
  VLMCKmer kmer(actual_kmer_size, 0, {});
  VLMCKmer prev_kmer(actual_kmer_size, 0, {});

  int alphabet_size = 4;
  // actual_kmer_size + 2 as 0 is used for the root, and end is used for
  // next-counts of the longest kmers.
  KMersPerLevel<uint64> counters{alphabet_size + 1, kmer_size + 1};

  uint64 counter;

  auto output_func = [&](const VLMCKmer &prev_kmer, int diff_pos, uint64 sum,
                         std::vector<uint64> &next_counts) {
    output_node(prev_kmer, diff_pos, sum, next_counts, container);
  };

  while (kmer_database.ReadNextKmer(kmer_api, counter)) {
    kmer_api.construct_vlmc_kmer(kmer);
    kmer.count = counter;

    if (prev_kmer.length != 0) {
      process_kmer(kmer, prev_kmer, counters, output_func, include_node, 1);
    }

    // The counts should be increased by this kmer's count
    increase_counts(counters, kmer);

    prev_kmer = kmer;
  }

  VLMCKmer epsilon(0, 0, {});
  process_kmer(epsilon, prev_kmer, counters, output_func, include_node, 1);
  output_root(counters[0], container);

  kmer_database.Close();
}

} // namespace vlmc
