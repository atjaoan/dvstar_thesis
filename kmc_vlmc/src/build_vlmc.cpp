#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>

#include <kmc_file.h>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

#include "kmer.hpp"
#include "similarity_pruning.hpp"
#include "support_pruning.hpp"

struct cli_arguments {
  std::string fasta_path;
  int min_count = 10;
  int max_depth = 15;
  double threshold = 3.9075;
};

std::string run_kmc(std::string fasta_path, const int kmer_size) {
  std::string kmc_db_name{"res"};
  std::string kmc_tmp{"tmp"};

  std::ostringstream stringStream;
  stringStream << "./kmc -b -ci1 -cs4294967295 ";
  stringStream << "-k" << kmer_size << " ";
  stringStream << "-fm " << fasta_path << " ";
  stringStream << kmc_db_name << " " << kmc_tmp;
  std::string command = stringStream.str();

  system(command.c_str());

  return kmc_db_name;
}

int main(int argc, char *argv[]) {
  CLI::App app{"Variable-length Markov chain construction construction using "
               "k-mer counter."};

  cli_arguments arguments{};

  app.add_option("-p,--fasta-path", arguments.fasta_path, "Path to fasta file.")
      ->required();
  app.add_option("-c,--min-count", arguments.min_count,
                 "Minimum count required for every k-mer in the tree.");
  app.add_option("-k,--threshold", arguments.threshold,
                 "Kullback-Leibler threshold.");
  app.add_option("-d,--max-depth", arguments.max_depth,
                 "Maximum depth for included k-mers.");

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  auto start = std::chrono::steady_clock::now();

  const int kmer_size = arguments.max_depth + 1;

  auto kmc_db_name = run_kmc(arguments.fasta_path, kmer_size);

  auto kmc_done = std::chrono::steady_clock::now();

  CKMCFile kmer_database;
  auto status = kmer_database.OpenForListing(kmc_db_name);

  if (!status) {
    std::cout << "opening file not successful" << std::endl;
    return EXIT_FAILURE;
  }

  // second parameter is RAM to use, should be parameter.
  kmer_sorter<31> sorter{ReverseKMerComparator<31>(), 128 * 1024 * 1024};

  auto include_node = [&](int length, size_t count) -> bool {
    return length <= arguments.max_depth && count >= arguments.min_count;
  };

  auto support_pruning_start = std::chrono::steady_clock::now();
  support_pruning(kmc_db_name, sorter, kmer_size, include_node);
  auto support_pruning_done = std::chrono::steady_clock::now();

  kmer_database.Close();

  auto sorting_start = std::chrono::steady_clock::now();
  sorter.sort();
  auto sorting_done = std::chrono::steady_clock::now();

  std::filesystem::path path{"out.txt"};
  std::ofstream ofs(path);

  auto keep_node = [&](double delta) -> bool {
    return delta <= arguments.threshold;
  };

  auto similarity_pruning_start = std::chrono::steady_clock::now();
  similarity_pruning(sorter, ofs, keep_node);
  auto similarity_pruning_done = std::chrono::steady_clock::now();

  ofs.close();

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