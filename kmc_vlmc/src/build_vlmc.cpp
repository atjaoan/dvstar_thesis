

#include <kmc_file.h>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

#include "build_vlmc.hpp"

int main(int argc, char *argv[]) {
  CLI::App app{"Variable-length Markov chain construction construction using "
               "k-mer counter."};

  cli_arguments arguments{};
  add_options(app, arguments);

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  std::filesystem::create_directories(arguments.tmp_path);
  configure_stxxl(arguments.tmp_path);

  if (arguments.mode == "build") {
    return build(arguments.fasta_path, arguments.max_depth, arguments.min_count,
                 arguments.threshold, arguments.out_path, arguments.tmp_path,
                 arguments.in_or_out_of_core);

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
  } else if (arguments.mode == "score") {
    negative_log_likelihood(arguments.fasta_path, arguments.tmp_path,
                            arguments.in_path, arguments.in_or_out_of_core,
                            arguments.max_depth);
  }

  std::filesystem::remove_all(arguments.tmp_path);

  return EXIT_SUCCESS;
}
