#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <kmc_file.h>

#include "kmer.hpp"
#include "support_pruning.hpp"

int main(int argc, char *argv[]) {
  uint32 min_count_to_set = 0;
  uint32 max_count_to_set = 0;
  const int kmer_size = 15;
  std::string input_file_name{"res"};


  CKMCFile kmer_database;
  auto status = kmer_database.OpenForListing(input_file_name);

  if (!status) {
    std::cout << "opening file not successful" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "test" << std::endl;

  // second parameter is RAM to use, should be parameter.
  kmer_sorter<kmer_size> sorter{KMerComparator<kmer_size>(), 64 * 1024 * 1024};
  support_pruning(kmer_database, sorter);

  sorter.sort();

  std::filesystem::path path{"out.txt"};
  std::ofstream ofs(path);

  while (!sorter.empty())
  {
    VLMCKmer kmer = *sorter;
    std::cout << kmer.to_string() << " " << kmer.count << " ";
    for (auto &c : kmer.next_symbol_counts) {
      std::cout << c << " ";
    }
    std::cout << std::endl;

    ++sorter;
  }

  ofs.close();
  kmer_database.Close();

  return EXIT_SUCCESS;
}