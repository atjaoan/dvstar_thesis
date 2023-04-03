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
#include "vlmc_from_kmers/kmer.hpp"

using VLMC_vector = container::VLMC_vector;
using VLMC_Indexing = container::VLMC_Indexing;

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

  std::function<double(VLMC_vector &, VLMC_vector &)> dist_func = [&](auto &left, auto &right) {
      return distance::dvstar<VLMC_vector>(left, right, background_order);
  };

  std::function<double(VLMC_Indexing &, VLMC_Indexing &)> dist_func_idx = [&](auto &left, auto &right) {
      return distance::dvstar<VLMC_Indexing>(left, right, background_order);
  };

  double error_tolerance = 1E-5;
};

TEST_F(DvstarTests, BackgroundOrderTest) {
  std::string test{"AAT"};
  vlmc::VLMCTranslator kmer{static_cast<int>(test.size())};
  if (!test.empty()) {
    kmer.from_string(test);
  }
  Kmer context = kmer.construct_vlmc_kmer();

  EXPECT_EQ("", distance::get_background_context(context.to_string(), 0));
  EXPECT_EQ("T", distance::get_background_context(context.to_string(), 1));
  EXPECT_EQ("AT", distance::get_background_context(context.to_string(), 2));
}

// Vector Tests
TEST_F(DvstarTests, Identity_vector) {
  VLMC_vector first_vlmc{first_bintree};

  double dist_vector = dist_func(first_vlmc, first_vlmc);
  EXPECT_NEAR(0.0, dist_vector, error_tolerance);
}

TEST_F(DvstarTests, EqualDistance_vector) {
  VLMC_vector first_vlmc{first_bintree};
  VLMC_vector second_vlmc{second_bintree};

  double dist_vector = dist_func(first_vlmc, second_vlmc);
  double old_dvstar_implementation = vlmc::dvstar(first_bintree, second_bintree, background_order);

  EXPECT_NEAR(old_dvstar_implementation, dist_vector, error_tolerance);
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
TEST_F(DvstarTests, Identity_indexing) {
  VLMC_Indexing first_vlmc{first_bintree};

  double dist_multi_vector = dist_func_idx(first_vlmc, first_vlmc);
  EXPECT_NEAR(0.0, dist_multi_vector, error_tolerance);
}

TEST_F(DvstarTests, Symmetry_indexing) {
  VLMC_Indexing first_vlmc{first_bintree};
  VLMC_Indexing second_vlmc{second_bintree};

  double dist_one_multi_vector = dist_func_idx(first_vlmc, second_vlmc);
  double dist_two_multi_vector = dist_func_idx(second_vlmc, first_vlmc);

  EXPECT_DOUBLE_EQ(dist_one_multi_vector, dist_two_multi_vector);
}

TEST_F(DvstarTests, multiple_runs_indexing) {
  VLMC_Indexing first_vlmc{first_bintree};
  VLMC_Indexing second_vlmc{second_bintree};
  size_t runs = 10; 
  double prev = dist_func_idx(first_vlmc, second_vlmc);
  for (int i = 0; i < runs; i++){
    double current = dist_func_idx(first_vlmc, second_vlmc);
    EXPECT_DOUBLE_EQ(prev, current); 
    prev = current; 
  }
}

TEST_F(DvstarTests, EqualDistance_indexing) {
  VLMC_Indexing first_vlmc{first_bintree};
  VLMC_Indexing second_vlmc{second_bintree};

  VLMC_vector first_vlmc_v{first_bintree};
  VLMC_vector second_vlmc_v{second_bintree};

  double dist_vector = dist_func(first_vlmc_v, second_vlmc_v);
  double dist_multi_vector = dist_func_idx(first_vlmc, second_vlmc);
  double old_dvstar_implementation = vlmc::dvstar(first_bintree, second_bintree, background_order);

  EXPECT_NEAR(dist_multi_vector, dist_vector, error_tolerance);
  EXPECT_NEAR(old_dvstar_implementation, dist_vector, error_tolerance);
  EXPECT_NEAR(old_dvstar_implementation, dist_multi_vector, error_tolerance);
}

TEST_F(DvstarTests, TestBackgroundOrder) {
  for (int order = 0; order < 5; order++){
    container::VLMC_vector first_vlmc_vector{first_bintree, order};
    container::VLMC_vector second_vlmc_vector{second_bintree, order};
    container::VLMC_Indexing first_vlmc_indexing{first_bintree, order};
    container::VLMC_Indexing second_vlmc_indexing{second_bintree, order};

    auto dist_vector = distance::dvstar<VLMC_vector>(first_vlmc_vector, second_vlmc_vector, order);
    auto dist_indexing = distance::dvstar<VLMC_Indexing>(first_vlmc_indexing, second_vlmc_indexing, order);
    auto old_dvstar_implementation = vlmc::dvstar(first_bintree, second_bintree, order);

    EXPECT_NEAR(old_dvstar_implementation, dist_vector, error_tolerance);
    EXPECT_NEAR(dist_vector, dist_indexing, error_tolerance);
  }
}
