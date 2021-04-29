#pragma once

#include <filesystem>
#include <random>
#include <string>

#include "kmer_container.hpp"

struct cli_arguments {
  std::string mode{"build"};
  std::filesystem::path in_path;
  std::filesystem::path tmp_path{};
  std::filesystem::path out_path{};
  int min_count = 10;
  int max_depth = 15;
  double threshold = 3.9075;
  std::string in_or_out_of_core{"internal"};
};

std::random_device rd;
std::mt19937 gen = std::mt19937{rd()};

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

std::unique_ptr<KmerContainer>
parse_kmer_container(const std::string &in_or_out_of_core) {
  std::cout << in_or_out_of_core << std::endl;
  if (in_or_out_of_core == "external") {
    return std::make_unique<OutOfCoreKmerContainer>();
  } else if (in_or_out_of_core == "internal") {
    return std::make_unique<InCoreKmerContainer>();
  } else {
    throw(std::invalid_argument(
        "parameter --out-or-in-of-core not 'internal' or 'external'"));
  }
}