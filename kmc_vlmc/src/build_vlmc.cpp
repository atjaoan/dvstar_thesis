#include <chrono>
#include <execution>
#include <filesystem>
#include <iostream>
#include <string>

#include <kmc_file.h>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>

#include "cli_helper.hpp"
#include "kmer_container.hpp"
#include "similarity_pruning.hpp"
#include "support_pruning.hpp"

std::filesystem::path run_kmc(std::filesystem::path &in_path,
                              const int kmer_size,
                              const std::filesystem::path &tmp_path) {
  std::string random_name = get_random_name("kmc");

  std::filesystem::path kmc_db_name = tmp_path / (random_name + "_res");
  std::filesystem::path kmc_tmp = tmp_path / (random_name + "_tmp");

  std::ostringstream stringStream;
  stringStream << "./kmc -b -ci1 -cs4294967295 ";
  stringStream << "-k" << kmer_size << " ";
  stringStream << "-fm " << in_path << " ";
  stringStream << kmc_db_name << " " << kmc_tmp;
  std::string command = stringStream.str();

  system(command.c_str());

  return kmc_db_name;
}

int main(int argc, char *argv[]) {
  CLI::App app{"Variable-length Markov chain construction construction using "
               "k-mer counter."};

  cli_arguments arguments{};

  app.add_option("-m,--mode", arguments.mode,
                 "Program mode, 'build', 'dump', or 'score'.");

  app.add_option("-p,--in-path", arguments.in_path, "Path to fasta file.")
      ->required();
  app.add_option("-t,--temp-path", arguments.tmp_path,
                 "Path to temporary folder for intermediate outputs.");
  app.add_option("-o,--out-path", arguments.out_path,
                 "Path to output file.  If no file is given, writes result to "
                 "standard out.");

  app.add_option("-c,--min-count", arguments.min_count,
                 "Minimum count required for every k-mer in the tree.");
  app.add_option("-k,--threshold", arguments.threshold,
                 "Kullback-Leibler threshold.");
  app.add_option("-d,--max-depth", arguments.max_depth,
                 "Maximum depth for included k-mers.");

  app.add_option(
      "-i, --in-or-out-of-core", arguments.in_or_out_of_core,
      "Specify 'internal' for in-core or 'external for out-of-core memory "
      "model.  Out of core is slower, but is not memory bound. ");

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  if (arguments.mode == "build") {
    std::filesystem::create_directories(arguments.tmp_path);

    auto start = std::chrono::steady_clock::now();

    const int kmer_size = arguments.max_depth + 1;

    auto kmc_db_name =
        run_kmc(arguments.in_path, kmer_size, arguments.tmp_path);

    auto kmc_done = std::chrono::steady_clock::now();

    CKMCFile kmer_database;
    auto status = kmer_database.OpenForListing(kmc_db_name);

    if (!status) {
      std::cout << "opening file not successful" << std::endl;
      return EXIT_FAILURE;
    }

    auto include_node = [&](int length, size_t count) -> bool {
      return length <= arguments.max_depth && count >= arguments.min_count;
    };

    auto container = parse_kmer_container(arguments.in_or_out_of_core);

    auto support_pruning_start = std::chrono::steady_clock::now();
    support_pruning<31>(kmer_database, *container, kmer_size, include_node);
    auto support_pruning_done = std::chrono::steady_clock::now();

    kmer_database.Close();

    auto sorting_start = std::chrono::steady_clock::now();
    (*container).sort();
    auto sorting_done = std::chrono::steady_clock::now();

    std::filesystem::path path =
        std::filesystem::path(arguments.in_path).stem();
    path += ".txt";

    std::ofstream file_stream(arguments.out_path, std::ios::binary);
    cereal::BinaryOutputArchive oarchive(file_stream);

    auto keep_node = [&](double delta) -> bool {
      return delta <= arguments.threshold;
    };

    auto similarity_pruning_start = std::chrono::steady_clock::now();
    similarity_pruning<31>(*container, oarchive, keep_node);
    auto similarity_pruning_done = std::chrono::steady_clock::now();

    if (!arguments.out_path.empty()) {
      file_stream.close();
    }

    std::chrono::duration<double> total_seconds =
        similarity_pruning_done - start;
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

    //    std::filesystem::remove_all(arguments.tmp_path);

  } else if (arguments.mode == "dump") {
    std::ifstream file_stream(arguments.in_path, std::ios::binary);
    cereal::BinaryInputArchive iarchive(file_stream);

    std::ostream *ofs = &std::cout;
    std::ofstream out_stream(arguments.out_path);

    if (arguments.out_path.empty()) {
      ofs = &std::cout;
    } else {
      ofs = &out_stream;
    }

    VLMCKmer kmer{};

    bool data_left_to_read;
    while (file_stream.peek() != EOF) {
      try {
        iarchive(kmer);
        kmer.output(*ofs);
      } catch (const cereal::Exception &e) {
        std::cout << (file_stream.peek() == EOF) << std::endl;
        std::cout << e.what() << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  return EXIT_SUCCESS;
}
