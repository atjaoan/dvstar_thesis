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

using VLMC_vector = container::VLMC_vector;
using Index_by_value = container::Index_by_value;

class DvstarTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_fasta{"NC_028367.1.fa"};
  std::filesystem::path second_fasta{"NC_045512.2.fa"};
  std::filesystem::path third_fasta{"NC_001497.2.fa"};

  std::filesystem::path first_bintree{"../data/test_VLMCs/sequences_1.bintree"};
  std::filesystem::path second_bintree{"../data/test_VLMCs/sequences_2.bintree"};
  std::filesystem::path third_bintree{"../data/test_VLMCs/sequences_3.bintree"};

  int background_order = 0;

  std::function<double(container::VLMC_Container &, container::VLMC_Container &)> dist_func = [&](auto &left, auto &right) {
      return distance::dvstar(left, right, background_order);
  };
};

// Vector Tests
TEST_F(DvstarTests, Identity_vector) {
  VLMC_vector first_vlmc{first_bintree};

  double dist_vector = dist_func(first_vlmc, first_vlmc);
  EXPECT_DOUBLE_EQ(0.0, dist_vector);
  double old_dvstar_implementation = vlmc::dvstar(first_bintree, first_bintree, background_order);
  EXPECT_DOUBLE_EQ(0.0, old_dvstar_implementation);
}

TEST_F(DvstarTests, EqualDistance_vector) {
  VLMC_vector first_vlmc{first_bintree};
  VLMC_vector second_vlmc{second_bintree};

  double dist_vector = dist_func(first_vlmc, second_vlmc);
  double old_dvstar_implementation = vlmc::dvstar(first_bintree, second_bintree, background_order);

  EXPECT_DOUBLE_EQ(old_dvstar_implementation, dist_vector);
}

TEST_F(DvstarTests, Symmetry_vector) {
  VLMC_vector first_vlmc{first_bintree};
  VLMC_vector second_vlmc{second_bintree};

  double dist_one_vector = dist_func(first_vlmc, second_vlmc);
  double dist_two_vector = dist_func(second_vlmc, first_vlmc);

  EXPECT_DOUBLE_EQ(dist_one_vector, dist_two_vector);
}

TEST_F(DvstarTests, multiple_runs_vector) {
  VLMC_vector first_vlmc{first_bintree};
  VLMC_vector second_vlmc{second_bintree};
  size_t runs = 10; 
  double prev = dist_func(first_vlmc, second_vlmc);
  for (int i = 0; i < runs; i++){
    double current = dist_func(first_vlmc, second_vlmc);
    EXPECT_DOUBLE_EQ(prev, current); 
    prev = current; 
  }
}

// Multi Vector Tests
TEST_F(DvstarTests, Identity_multivec) {
  Index_by_value first_vlmc{first_bintree};

  double dist_multi_vector = dist_func(first_vlmc, first_vlmc);
  EXPECT_DOUBLE_EQ(0.0, dist_multi_vector);
  double old_dvstar_implementation = vlmc::dvstar(first_bintree, first_bintree, background_order);
  EXPECT_DOUBLE_EQ(0.0, old_dvstar_implementation);
}

TEST_F(DvstarTests, Symmetry_multivec) {
  Index_by_value first_vlmc{first_bintree};
  Index_by_value second_vlmc{second_bintree};

  double dist_one_multi_vector = dist_func(first_vlmc, second_vlmc);
  double dist_two_multi_vector = dist_func(second_vlmc, first_vlmc);

  EXPECT_DOUBLE_EQ(dist_one_multi_vector, dist_two_multi_vector);
}

TEST_F(DvstarTests, multiple_runs_multivec) {
  Index_by_value first_vlmc{first_bintree};
  Index_by_value second_vlmc{second_bintree};
  size_t runs = 10; 
  double prev = dist_func(first_vlmc, second_vlmc);
  for (int i = 0; i < runs; i++){
    double current = dist_func(first_vlmc, second_vlmc);
    EXPECT_DOUBLE_EQ(prev, current); 
    prev = current; 
  }
}

TEST_F(DvstarTests, EqualDistance_multivec) {
  Index_by_value first_vlmc{first_bintree};
  Index_by_value second_vlmc{second_bintree};

  VLMC_vector first_vlmc_v{first_bintree};
  VLMC_vector second_vlmc_v{second_bintree};

  double dist_vector = dist_func(first_vlmc_v, second_vlmc_v);
  double dist_multi_vector = dist_func(first_vlmc, second_vlmc);
  double old_dvstar_implementation = vlmc::dvstar(first_bintree, second_bintree, background_order);

  EXPECT_DOUBLE_EQ(dist_multi_vector, dist_vector);
  EXPECT_DOUBLE_EQ(old_dvstar_implementation, dist_vector);
  EXPECT_DOUBLE_EQ(old_dvstar_implementation, dist_multi_vector);
}
