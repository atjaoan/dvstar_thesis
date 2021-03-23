#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <kmc_file.h>

#include "support_pruning.hpp"

int main(int argc, char *argv[]) {
  uint32 min_count_to_set = 0;
  uint32 max_count_to_set = 0;
  size_t kmer_size = 15;
  std::string input_file_name{"res"};

  std::filesystem::path path{"out.txt"};
  std::ofstream ofs(path);

  CKMCFile kmer_database;
  auto status = kmer_database.OpenForListing(input_file_name);

  if (!status) {
    std::cout << "opening file not successful" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "test" << std::endl;


  ofs.close();
  kmer_database.Close();

  return EXIT_SUCCESS;
}