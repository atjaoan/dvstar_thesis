#pragma once

#include <filesystem>
#include <string>

#include "cli_helper.hpp"

std::filesystem::path run_kmc(const std::filesystem::path &fasta_path,
                              const int kmer_size,
                              const std::filesystem::path &tmp_path,
                              const int exclude_infrequent_count = 2) {
  std::string random_name = get_random_name("kmc_");

  std::filesystem::path kmc_db_name = tmp_path / (random_name + "_res");
  std::filesystem::path kmc_tmp = tmp_path / (random_name + "_tmp");

  std::ostringstream stringStream;
  stringStream << "./kmc -b -ci" << exclude_infrequent_count
               << " -cs4294967295 ";
  stringStream << "-k" << kmer_size << " ";
  stringStream << "-fm " << fasta_path << " ";
  stringStream << kmc_db_name << " " << kmc_tmp;
  std::string command = stringStream.str();

  system(command.c_str());

  return kmc_db_name;
}
