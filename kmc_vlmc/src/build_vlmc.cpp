#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <kmc_file.h>

#include "kmer.hpp"
#include "support_pruning.hpp"

int main(int argc, char *argv[]) {
  const int kmer_size = 10;
  std::string input_file_name{"res"};


  CKMCFile kmer_database;
  auto status = kmer_database.OpenForListing(input_file_name);

  if (!status) {
    std::cout << "opening file not successful" << std::endl;
    return EXIT_FAILURE;
  }

  // second parameter is RAM to use, should be parameter.
  kmer_sorter<kmer_size> sorter{ReverseKMerComparator<kmer_size>(), 64 * 1024 * 1024};
  support_pruning(kmer_database, sorter);

  sorter.sort();

  std::filesystem::path path{"out.txt"};
  std::ofstream ofs(path);

  while (!sorter.empty())
  {
    VLMCKmer kmer = *sorter;
    kmer.output(ofs);

    ++sorter;
  }

  ofs.close();
  kmer_database.Close();

  return EXIT_SUCCESS;
}