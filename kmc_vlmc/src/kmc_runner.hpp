#pragma once

#include <filesystem>
#include <string>

#include "cli_helper.hpp"

std::filesystem::path run_kmc(const std::filesystem::path &fasta_path,
                              const int kmer_size,
                              const std::filesystem::path &tmp_path,
                              const std::string &in_or_out_of_core,
                              const int exclude_infrequent_count = 2) {
  std::string random_name = get_random_name("kmc_");

  std::filesystem::path kmc_db_name = tmp_path / (random_name + "_res");
  std::filesystem::path kmc_db_sorted_name =
      tmp_path / (random_name + "_res_sorted");
  std::filesystem::path kmc_db_sorted_name_with_suffix =
      tmp_path / (random_name + "_res_sorted.kmc_suf");
  std::filesystem::path kmc_tmp = tmp_path / (random_name + "_tmp");

  std::filesystem::create_directories(kmc_tmp);

  auto kmc_run = std::chrono::steady_clock::now();

  std::ostringstream kmc_run_stream;
  kmc_run_stream << "./kmc -b -ci" << exclude_infrequent_count
                 << " -cs4294967295 ";
  if (in_or_out_of_core == "internal") {
    kmc_run_stream << "-r ";
  }
  kmc_run_stream << "-k" << kmer_size << " ";
  kmc_run_stream << "-fm " << fasta_path << " ";
  kmc_run_stream << kmc_db_name << " " << kmc_tmp;
  std::string kmc_run_command = kmc_run_stream.str();

  system(kmc_run_command.c_str());
  auto kmc_run_done = std::chrono::steady_clock::now();

  auto kmc_sort = std::chrono::steady_clock::now();

  std::ostringstream kmc_sort_stream;
  kmc_sort_stream << "./kmc_tools transform ";
  kmc_sort_stream << kmc_db_name;
  kmc_sort_stream << " sort ";
  kmc_sort_stream << kmc_db_sorted_name;
  std::string kmc_sort_command = kmc_sort_stream.str();

  system(kmc_sort_command.c_str());

  auto kmc_sort_done = std::chrono::steady_clock::now();

  std::chrono::duration<double> kmc_run_seconds = kmc_run_done - kmc_run;
  std::cout << "KMC count time: " << kmc_run_seconds.count() << "s\n";

  std::chrono::duration<double> kmc_sort_seconds = kmc_sort_done - kmc_sort;
  std::cout << "KMC sort time: " << kmc_sort_seconds.count() << "s\n";

  if (std::filesystem::exists(kmc_db_sorted_name_with_suffix)) {
    return kmc_db_sorted_name;
  } else {
    return kmc_db_name;
  }
}
