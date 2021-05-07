#pragma once

#include <filesystem>

#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>
#include <chrono>
#include <execution>
#include <filesystem>
#include <iostream>
#include <string>

#include "cli_helper.hpp"
#include "kmc_runner.hpp"
#include "kmer_container.hpp"
#include "negative_log_likelihood.hpp"

#include "similarity_pruning.hpp"
#include "support_pruning.hpp"

void configure_stxxl(const std::filesystem::path &tmp_path) {
  std::string random_name = get_random_name("stxxl_");
  // get uninitialized config singleton
  stxxl::config *cfg = stxxl::config::get_instance();
  std::filesystem::path stxxl_disk_path = tmp_path / random_name;
  // create a disk_config structure.
  // allocate 6GB
  stxxl::disk_config disk1(stxxl_disk_path, uint64(6000) * 1024 * 1024,
                           "linuxaio unlink");

  cfg->add_disk(disk1);
}

int build(const std::filesystem::path &fasta_path, const int max_depth,
          const int min_count, const double threshold,
          const std::filesystem::path &out_path,
          const std::filesystem::path &tmp_path,
          const std::string &in_or_out_of_core) {
  auto start = std::chrono::steady_clock::now();

  const int kmer_size = max_depth + 1;

  auto kmc_db_name =
      run_kmc(fasta_path, kmer_size, tmp_path, in_or_out_of_core, 2);

  auto kmc_done = std::chrono::steady_clock::now();

  CKMCFile kmer_database;
  auto status = kmer_database.OpenForListing(kmc_db_name);

  if (!status) {
    std::cout << "opening file not successful" << std::endl;
    return EXIT_FAILURE;
  }

  auto include_node = [&](int length, size_t count) -> bool {
    return length <= max_depth && count >= min_count;
  };

  auto container =
      parse_kmer_container<ReverseKMerComparator<31>>(in_or_out_of_core);

  auto support_pruning_start = std::chrono::steady_clock::now();
  support_pruning<31>(kmer_database, container, kmer_size, include_node, in_or_out_of_core);
  auto support_pruning_done = std::chrono::steady_clock::now();

  kmer_database.Close();

  auto sorting_start = std::chrono::steady_clock::now();

  container->sort();
  auto sorting_done = std::chrono::steady_clock::now();

  std::filesystem::path path = std::filesystem::path(fasta_path).stem();
  path += ".txt";

  std::ofstream file_stream(out_path, std::ios::binary);
  auto similarity_pruning_start = std::chrono::steady_clock::now();
  {
    cereal::BinaryOutputArchive oarchive(file_stream);

    auto keep_node = [&](double delta) -> bool { return delta <= threshold; };

    similarity_pruning<31>(container, oarchive, keep_node);
  }
  auto similarity_pruning_done = std::chrono::steady_clock::now();

  if (!out_path.empty()) {
    file_stream.close();
  }

  std::chrono::duration<double> total_seconds = similarity_pruning_done - start;
  std::cout << "Total time: " << total_seconds.count() << "s\n";

  std::chrono::duration<double> kmc_seconds = kmc_done - start;
  std::cout << "KMC time: " << kmc_seconds.count() << "s\n";

  std::chrono::duration<double> support_seconds =
      support_pruning_done - support_pruning_start;
  std::cout << "Support pruning time: " << support_seconds.count() << "s\n";

  std::chrono::duration<double> sort_seconds = sorting_done - sorting_start;
  std::cout << "Sorting time: " << sort_seconds.count() << "s\n";

  std::chrono::duration<double> similarity_seconds =
      similarity_pruning_done - sorting_done;
  std::cout << "Similarity pruning time: " << similarity_seconds.count()
            << "s\n";

  return EXIT_SUCCESS;
}