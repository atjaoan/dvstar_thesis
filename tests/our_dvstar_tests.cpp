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
};


TEST_F(OurDvstarTests, Identity) {
  vlmc_c first_vlmc{first_bintree};

  double dist = distance::dvstar(first_vlmc, first_vlmc);

  EXPECT_FLOAT_EQ(dist, 0.0);
  // double old_dist = vlmc::dvstar(first_bintree, first_bintree, 0);
  // EXPECT_FLOAT_EQ(dist, old_dist);
}