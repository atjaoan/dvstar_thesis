#pragma once

#include <filesystem>
#include <memory>
#include <random>
#include <string>

#include "kmer_container.hpp"

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

namespace vlmc {

enum Mode {
  build,
  score_sequence,
  dump,
  bic,
  build_from_kmc_db,
  dvstar_dissimliarity
};

struct cli_arguments {
  Mode mode{Mode::build};
  std::filesystem::path fasta_path;
  std::filesystem::path in_path;
  std::filesystem::path to_path;
  std::filesystem::path tmp_path{"./tmp"};
  std::filesystem::path out_path{};
  int min_count = 10;
  int max_depth = 15;
  double threshold = 3.9075;
  double pseudo_count_amount = 1.0;
  Core in_or_out_of_core{Core::in};
};
void add_options(CLI::App &app, cli_arguments &arguments) {
  std::map<std::string, Mode> mode_map{
      {"build", Mode::build},
      {"score", Mode::score_sequence},
      {"dump", Mode::dump},
      {"bic", Mode::bic},
      {"build-from-kmc-db", Mode::build_from_kmc_db},
      {"dvstar", Mode::dvstar_dissimliarity}};

  std::map<std::string, Core> core_map{
      {"internal", Core::in}, {"external", Core::out}, {"hash", Core::hash}};

  app.add_option(
         "-m,--mode", arguments.mode,
         "Program mode, 'build', 'build-from-kmc-db', 'dump', 'score', or "
         "'dvstar'.  For "
         "build-from-kmc-db, the kmc db needs to include all k-mers (not in "
         "canonical form), with no minimum count cutoff.  The length of the "
         "k-mers needs to be set to 1 more than the maximum depth of the VLMC. "
         " The kmc db also has to be sorted.")
      ->transform(CLI::CheckedTransformer(mode_map, CLI::ignore_case));

  app.add_option(
      "-p,--fasta-path", arguments.fasta_path,
      "Path to fasta file.  Required for 'build' and 'score' modes.");

  app.add_option(
      "--in-path", arguments.in_path,
      "Path to saved tree file or kmc db file.  Required for "
      "'build-from-kmc-db', 'dump', 'score', and 'dvstar' modes.  For "
      "'build-from-kmc-db', the kmc db file needs to be supplied "
      "without the file extension.");

  app.add_option("--to-path", arguments.to_path,
                 "Path to saved tree file.  Required for 'dvstar' mode.");

  app.add_option("-o,--out-path", arguments.out_path,
                 "Path to output file.  The VLMCs are stored as binary, and "
                 "can be read by the 'dump' or 'score' modes.  Required for "
                 "'build' and 'dump' modes.");

  app.add_option("-t,--temp-path", arguments.tmp_path,
                 "Path to temporary folder for the external memory algorithms. "
                 " For good performance, this needs to be on a local machine.  "
                 "For sorting, at least 2GB will be allocated to this path.  "
                 "Defaults to ./tmp");

  app.add_option("-c,--min-count", arguments.min_count,
                 "Minimum count required for every k-mer in the tree.");
  app.add_option("-k,--threshold", arguments.threshold,
                 "Kullback-Leibler threshold.");

  app.add_option("-d,--max-depth", arguments.max_depth,
                 "Maximum depth for included k-mers.");

  app.add_option("-a,--pseudo-count-amount", arguments.pseudo_count_amount,
                 "Size of pseudo count for probability estimation. See e.g. "
                 "https://en.wikipedia.org/wiki/Additive_smoothing .");

  app.add_option(
         "-i, --in-or-out-of-core", arguments.in_or_out_of_core,
         "Specify 'internal' for in-core or 'external for out-of-core memory "
         "model.  Out of core is slower, but is not memory bound. ")
      ->transform(CLI::CheckedTransformer(core_map, CLI::ignore_case));
}

static std::random_device rd;

static std::mt19937 gen = std::mt19937{rd()};
std::string get_random_name(const std::string &start) {
  std::stringstream ss;
  ss << start;

  // Generate 20 characters of random lower-case letters
  std::uniform_int_distribution<> distrib(97, 122);
  for (size_t i = 0; i < 20; i++) {
    auto rand = distrib(gen);
    ss << char(rand);
  }

  return ss.str();
}

template <class Comparator>
std::shared_ptr<KmerContainer<Comparator>>
parse_kmer_container(const Core &in_or_out_of_core) {
  if (in_or_out_of_core == Core::out) {
    return std::make_shared<OutOfCoreKmerContainer<Comparator>>();
  } else if (in_or_out_of_core == Core::in) {
    return std::make_shared<InCoreKmerContainer<Comparator>>();
  } else if (in_or_out_of_core == Core::hash) {
    return std::make_shared<HashMapKmerContainer<Comparator>>();
  } else {
    throw(std::invalid_argument(
        "parameter --out-or-in-of-core not 'internal' or 'external'"));
  }
}

} // namespace vlmc
