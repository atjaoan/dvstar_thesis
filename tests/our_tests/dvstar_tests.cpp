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

class DvstarTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_fasta{"NC_028367.1.fa"};
  std::filesystem::path second_fasta{"NC_045512.2.fa"};
  std::filesystem::path third_fasta{"NC_001497.2.fa"};

  std::filesystem::path first_bintree{"../data/test_VLMCs/sequences_1.bintree"};
  std::filesystem::path second_bintree{"../data/test_VLMCs/sequences_2.bintree"};
  std::filesystem::path third_bintree{"../data/test_VLMCs/sequences_3.bintree"};
};


TEST_F(DvstarTests, Identity_vector) {
  container::VLMC_vector first_vlmc{first_bintree};

  double dist = distance::dvstar(first_vlmc, first_vlmc);
  EXPECT_DOUBLE_EQ(dist, 0.0);
  double dist_real = vlmc::dvstar(first_bintree, first_bintree, 0);
  EXPECT_DOUBLE_EQ(dist_real, 0.0);
}

TEST_F(DvstarTests, EqualDistance_vector) {
  container::VLMC_vector first_vlmc{first_bintree};
  container::VLMC_vector second_vlmc{second_bintree};

  double dist_our = distance::dvstar(first_vlmc, second_vlmc);
  double dist_real = vlmc::dvstar(first_bintree, second_bintree, 0);

  EXPECT_DOUBLE_EQ(dist_our, dist_real);
}

TEST_F(DvstarTests, Symmetry_vector) {
  container::VLMC_vector first_vlmc{first_bintree};
  container::VLMC_vector second_vlmc{second_bintree};

  double dist_one = distance::dvstar(first_vlmc, second_vlmc);
  double dist_two = distance::dvstar(second_vlmc, first_vlmc);

  EXPECT_DOUBLE_EQ(dist_one, dist_two);
}

TEST_F(DvstarTests, Identity_multivec) {
  container::VLMC_multi_vector first_vlmc{first_bintree};

  double dist = distance::dvstar(first_vlmc, first_vlmc);
  EXPECT_DOUBLE_EQ(dist, 0.0);
  double dist_real = vlmc::dvstar(first_bintree, first_bintree, 0);
  EXPECT_DOUBLE_EQ(dist_real, 0.0);
}

TEST_F(DvstarTests, EqualDistance_multivec) {
  container::VLMC_multi_vector first_vlmc{first_bintree};
  container::VLMC_multi_vector second_vlmc{second_bintree};

  container::VLMC_vector first_vlmc_v{first_bintree};
  container::VLMC_vector second_vlmc_v{second_bintree};

  double dist_prev = distance::dvstar(first_vlmc_v, second_vlmc_v);
  double dist_our = distance::dvstar(first_vlmc, second_vlmc);
  double dist_real = vlmc::dvstar(first_bintree, second_bintree, 0);

  EXPECT_DOUBLE_EQ(dist_our, dist_prev);
  EXPECT_DOUBLE_EQ(dist_our, dist_real);
}

TEST_F(DvstarTests, Symmetry_multivec) {
  container::VLMC_multi_vector first_vlmc{first_bintree};
  container::VLMC_multi_vector second_vlmc{second_bintree};

  double dist_one = distance::dvstar(first_vlmc, second_vlmc);
  double dist_two = distance::dvstar(second_vlmc, first_vlmc);

  EXPECT_DOUBLE_EQ(dist_one, dist_two);
}
