#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

#include <kmc_file.h>

#include "vlmc_from_kmers/build_vlmc.hpp"
#include "read_helper.hpp"

using namespace vlmc;

class MethodRegression : public ::testing::Test {
protected:
  void SetUp() override {}
};

std::unordered_map<std::string, VLMCKmer>
to_map(const std::filesystem::path &run_path) {
  std::unordered_map<std::string, VLMCKmer> map{};

  std::ifstream run_stream{run_path};

  std::ifstream file_stream(run_path, std::ios::binary);
  {
    cereal::BinaryInputArchive iarchive(file_stream);

    VLMCKmer kmer{};
    while (file_stream.peek() != EOF) {
      iarchive(kmer);
      auto kmer_string = kmer.to_string();
      map[kmer_string] = kmer;
    }
  }

  file_stream.close();

  return map;
}

void test(vlmc::Core one, vlmc::Core two, int max_depth = 3) {
  std::filesystem::path fasta_path{"../python-prototype/NC_022098.1.fasta"};

  std::filesystem::path tmp_path = std::filesystem::temp_directory_path();

  std::filesystem::create_directories(tmp_path);
  vlmc::configure_stxxl(tmp_path);

  std::filesystem::path run_one_path{"out-one.tree"};
  int exit_one_code = vlmc::build_vlmc(fasta_path, max_depth, 2, 3.9075,
                                       run_one_path, tmp_path, one);

  std::filesystem::path run_two_path{"out-two.tree"};
  int exit_two_code = vlmc::build_vlmc(fasta_path, max_depth, 2, 3.9075,
                                       run_two_path, tmp_path, two);

  auto one_map = to_map(run_one_path);
  auto two_map = to_map(run_two_path);

  for (auto &[ctx, two_v] : two_map) {
    ASSERT_NE(one_map.find(ctx), one_map.end())
        << ctx << " " << two_v.divergence;
  }

  for (auto &[ctx, one_v] : one_map) {
    if (ctx == "") {
      continue;
    }

    ASSERT_NE(two_map.find(ctx), two_map.end())
        << ctx << " " << one_v.divergence;

    if (two_map.find(ctx) == two_map.end()) {
      continue;
    }

    auto two_v = two_map[ctx];

    if (one_v.divergence != two_v.divergence && two_v.divergence != -1 &&
        one_v.divergence != -1) {
      one_v.output(std::cout);
      two_v.output(std::cout);
    }

    ASSERT_EQ(one_v, two_v);
    ASSERT_EQ(one_v.count, two_v.count) << ctx;
    if (two_v.divergence != -1 && one_v.divergence != -1) {
      ASSERT_EQ(one_v.divergence, two_v.divergence) << ctx;
    }
    for (int i = 0; i < 4; i++) {
      ASSERT_EQ(one_v.next_symbol_counts[i], two_v.next_symbol_counts[i])
          << ctx;
    }
  }
}

// TEST_F(MethodRegression, InAndOut) {
//   test(vlmc::Core::in, vlmc::Core::out);
// }

TEST_F(MethodRegression, InAndHash) {
  for (int i = 1; i < 25; i++) {
    test(vlmc::Core::in, vlmc::Core::hash, i);
  }
}
