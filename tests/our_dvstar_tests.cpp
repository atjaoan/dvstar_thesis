#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

//Our files 
#include "vlmc_container.hpp"
#include "distances/dvstar.hpp"
#include "vlmc_from_kmers/dvstar.hpp"

//Original implementation files
#include <kmc_file.h>
#include "vlmc_from_kmers/dvstar.hpp"
#include "vlmc_from_kmers/build_vlmc.hpp"

// using namespace dvstarnamespace;

using vlmc_c = container::VLMC_vector;

class OurDvstarTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_fasta{"NC_028367.1.fa"};
  std::filesystem::path second_fasta{"NC_045512.2.fa"};
  std::filesystem::path third_fasta{"NC_001497.2.fa"};

  std::filesystem::path first_bintree{"../data/test_VLMCs/sequences_1.bintree"};
  std::filesystem::path second_bintree{"../data/test_VLMCs/sequences_2.bintree"};
  std::filesystem::path third_bintree{"../data/test_VLMCs/sequences_3.bintree"};

  vlmc_c first_vlmc{first_bintree};
  vlmc_c second_vlmc{second_bintree};
};


TEST_F(OurDvstarTests, Identity) {
  std::filesystem::path tmp_path = std::filesystem::temp_directory_path();

  std::filesystem::create_directories(tmp_path);
  vlmc::configure_stxxl(tmp_path);

  std::filesystem::path run_one_path{"out-one.bintree"};
  int exit_one_code = vlmc::build_vlmc(first_fasta, 6, 2, 3.9075,
                                       run_one_path, tmp_path, vlmc::Core::in);

  std::filesystem::path run_two_path{"out-two.bintree"};
  int exit_two_code = vlmc::build_vlmc(first_fasta, 6, 2, 3.9075,
                                       run_two_path, tmp_path, vlmc::Core::in);

  double dist_test = vlmc::dvstar(run_one_path, run_two_path, 0);

  double dist = distance::dvstar(first_vlmc, first_vlmc);
  double dist_real = vlmc::dvstar(first_bintree, first_bintree, 0);
  EXPECT_DOUBLE_EQ(dist_test, 0.0);
}

TEST_F(OurDvstarTests, EqualDistance) {

  double dist_our = distance::dvstar(first_vlmc, second_vlmc);
  double dist_real = vlmc::dvstar(first_bintree, second_bintree, 0);

  EXPECT_DOUBLE_EQ(dist_our, dist_real);
}