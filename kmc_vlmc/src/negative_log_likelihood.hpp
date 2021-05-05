#pragma once

#include <filesystem>
#include <memory>

#include "kmc_runner.hpp"
#include "kmer_container.hpp"
#include "kmers_per_level.hpp"

void load_kmers(const std::filesystem::path &vlmc_path,
                KmerContainer<KMerComparator<31>> &container) {
  std::ifstream file_stream(vlmc_path, std::ios::binary);
  {
    cereal::BinaryInputArchive iarchive(file_stream);

    VLMCKmer kmer{};
    while (file_stream.peek() != EOF) {
      iarchive(kmer);
      container.push(kmer);
    }
  }

  file_stream.close();
}

std::tuple<double, bool, bool> compare_kmers(VLMCKmer &vlmc_kmer,
                                             VLMCKmer &kmc_kmer,
                                             const int prefix_length) {
  // Need to find the vlmc kmer that best matches the given kmer
  auto diff_pos = VLMCKmer::get_first_differing_position(vlmc_kmer, kmc_kmer, 0,
                                                         prefix_length);
  // If the two kmers differ anywhere (that isn't the last position),
  // they don't match.
  // If they don't differ anywhere in the vlmc kmer, we've found a match.
  if (diff_pos == -1) {
    auto char_idx = kmc_kmer.char_pos(kmc_kmer.length - 1);

    // Remember the pseudo-counts.
    double probability = double(vlmc_kmer.next_symbol_counts[char_idx] + 1) /
                         double(vlmc_kmer.count + 4);
    double log_likelihood = std::log(probability) * double(kmc_kmer.count);

    return {log_likelihood, true, false};
  } else {
    auto kmc_char_idx = kmc_kmer.char_pos(diff_pos + prefix_length);
    auto vlmc_char_idx = vlmc_kmer.char_pos(diff_pos);

    if (kmc_char_idx < vlmc_char_idx) {
      // Advance the kmc kmer, but save it as haven't found
      // find its matching vlmc kmer yet.
      return {0.0, true, true};
    } else {
      // Advance to the next vlmc kmer
      return {0.0, false, false};
    }
  }
}

std::tuple<double, stxxl::vector<VLMCKmer>>
score_kmers(KmerContainer<KMerComparator<31>> &container,
            const stxxl::vector<VLMCKmer> &kmers, const int kmer_length,
            const int prefix_length) {
  if (kmers.empty()) {
    return {0.0, {}};
  }
  stxxl::vector<VLMCKmer> next_level_kmers{};
  double log_likelihood = 0.0;

  bool next_please = false;
  bool kmc_done = false;
  int kmers_idx = 0;

  VLMCKmer kmc_kmer = kmers[kmers_idx];

  container.for_each([&](VLMCKmer &vlmc_kmer) {
    if (kmer_length != vlmc_kmer.length || kmc_done) {
      // Iterate the vlmc-kmers one depth at a time.
      return;
    }

    do {
      auto [log_likelihood_, next_please_, save_kmer] =
          compare_kmers(vlmc_kmer, kmc_kmer, prefix_length);
      if (save_kmer) {
        next_level_kmers.push_back(kmc_kmer);
      }
      next_please = next_please_;
      log_likelihood += log_likelihood_;

      if (next_please) {
        if (kmers.size() <= kmers_idx + 1) {
          kmc_done = true;
          return;
        } else {
          kmers_idx++;
          kmc_kmer = kmers[kmers_idx];
        }
      }
    } while (next_please);
  });

  return {log_likelihood, next_level_kmers};
}

std::tuple<double, stxxl::vector<VLMCKmer>, size_t>
score_kmers_full_length(KmerContainer<KMerComparator<31>> &container,
                        const std::string &kmc_db_name,
                        const int actual_kmer_size) {
  CKMCFile kmer_database;
  auto status = kmer_database.OpenForListing(kmc_db_name);

  if (!status) {
    std::cout << "opening file not successful" << std::endl;
    throw std::invalid_argument("internal: kmc db path does not open.");
  }

  VLMCTranslator kmer_api(actual_kmer_size + 1);
  uint64 counter;
  VLMCKmer kmc_kmer(actual_kmer_size + 1, 0, {});

  status = kmer_database.ReadNextKmer(kmer_api, counter);
  kmc_kmer = kmer_api.construct_vlmc_kmer();
  kmc_kmer.count = counter;

  if (!status) {
    return {0.0, {}, 0};
  }

  VLMCKmer prev_vlmc_kmer{};
  double log_likelihood = 0.0;
  size_t sequence_length = counter;

  stxxl::vector<VLMCKmer> next_level_kmers{};

  bool next_please = false;

  container.for_each([&](VLMCKmer &vlmc_kmer) {
    if (actual_kmer_size != vlmc_kmer.length) {
      // Iterate the vlmc-kmers one depth at a time.
      return;
    }

    do {
      auto [log_likelihood_, next_please_, save_kmer] =
          compare_kmers(vlmc_kmer, kmc_kmer, 0);

      if (save_kmer) {
        next_level_kmers.push_back(kmc_kmer);
      }

      next_please = next_please_;
      log_likelihood += log_likelihood_;

      if (next_please) {
        auto status = kmer_database.ReadNextKmer(kmer_api, counter);
        if (!status) {
          return;
        }
        kmc_kmer = kmer_api.construct_vlmc_kmer();
        kmc_kmer.count = counter;
        sequence_length += counter;
      }
    } while (next_please);
  });

  return {log_likelihood, next_level_kmers, sequence_length};
}

double score(KmerContainer<KMerComparator<31>> &container,
             const std::string &kmc_db_name, const int actual_kmer_size) {

  auto [full_log_likelihood, next_level_kmers, sequence_length] =
      score_kmers_full_length(container, kmc_db_name, actual_kmer_size);

  double log_likelihood = full_log_likelihood;

  int kmer_length = actual_kmer_size - 1;
  int prefix_length = 1;

  while (!next_level_kmers.empty()) {
    auto [log_likelihood_, next_level_kmers_] =
        score_kmers(container, next_level_kmers, kmer_length, prefix_length);

    next_level_kmers = next_level_kmers_;

    log_likelihood += log_likelihood_;

    kmer_length--;
    prefix_length++;
  }

  return -log_likelihood / double(sequence_length);
}

double negative_log_likelihood(const std::filesystem::path &fasta_path,
                               const std::filesystem::path &tmp_path,
                               const std::filesystem::path &vlmc_path,
                               const std::string &in_or_out_of_core,
                               const int actual_kmer_size) {
  auto kmer_container =
      parse_kmer_container<KMerComparator<31>>(in_or_out_of_core);
  load_kmers(vlmc_path, *kmer_container);

  kmer_container->sort();

  kmer_container->for_each([](VLMCKmer &kmer) { kmer.output(std::cout); });

  auto kmc_db_name =
      run_kmc(fasta_path, actual_kmer_size + 1, tmp_path, in_or_out_of_core, 1);

  double score_ = score(*kmer_container, kmc_db_name, actual_kmer_size);
  std::cout << "nll: " << score_ << std::endl;
  return score_;
}