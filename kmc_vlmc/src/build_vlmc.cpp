#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <kmc_file.h>

#include "kmer.hpp"
#include "similarity_pruning.hpp"
#include "support_pruning.hpp"

int main(int argc, char *argv[]) {
  const int kmer_size = 10;

  std::string kmc_db_name{"res"};
  std::string kmc_tmp{"tmp"};
  std::string fasta_path{"../python-prototype/NC_022098.1.fasta"};

  std::ostringstream stringStream;
  stringStream << "./kmc -b -ci1 -cs4294967295 ";
  stringStream << "-k" << kmer_size << " ";
  stringStream << "-fm " << fasta_path << " ";
  stringStream << kmc_db_name << " " << kmc_tmp;
  std::string command = stringStream.str();

  system(command.c_str());

  CKMCFile kmer_database;
  auto status = kmer_database.OpenForListing(kmc_db_name);

  if (!status) {
    std::cout << "opening file not successful" << std::endl;
    return EXIT_FAILURE;
  }


  // second parameter is RAM to use, should be parameter.
  kmer_sorter<kmer_size> sorter{ReverseKMerComparator<kmer_size>(),
                                64 * 1024 * 1024};
  support_pruning(kmer_database, sorter);
  kmer_database.Close();

  sorter.sort();

  std::filesystem::path path{"out.txt"};
  std::ofstream ofs(path);

  similarity_pruning(sorter, ofs);

  ofs.close();

  return EXIT_SUCCESS;
}