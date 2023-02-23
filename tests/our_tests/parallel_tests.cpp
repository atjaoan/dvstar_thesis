#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <cstdlib>
#include <filesystem>
#include <fstream>

//Our files 
#include "vlmc_container.hpp"
#include "distances/dvstar.hpp"
#include "vlmc_from_kmers/dvstar.hpp"
#include "optimize_cmp_kmers/get_cluster.hpp"
#include "optimize_cmp_kmers/calc_dists.hpp"

using matrix_t  = Eigen::MatrixXd;
using vlmc_t = container::VLMC_Container;

class ParallelTest : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_directory{"../data/test_VLMCs"};
  std::filesystem::path second_bintree{"../data/test_VLMCs"};

  size_t background_order = 0; 

  std::function<double(vlmc_t &, vlmc_t &)> distance_function = [&](auto &left, auto &right) {
      return distance::dvstar(left, right, background_order);
  };
};

TEST_F(ParallelTest, SequentialEqParallel) {
  // auto distance_function = distance::dvstar;
  auto cluster = cluster::get_cluster<container::VLMC_vector>(first_directory, 1);

  matrix_t distances_parallel{cluster.size(), cluster.size()};
  matrix_t distances_sequantial{cluster.size(), cluster.size()};

  auto fun_parallel = [&](size_t start_index, size_t stop_index) {
    calculate::calculate_full_slice<container::VLMC_vector>(start_index, stop_index, distances_parallel,
                           cluster, cluster, distance_function);
  };

  auto fun_sequential = [&](size_t start_index, size_t stop_index) {
    calculate::calculate_full_slice<container::VLMC_vector>(start_index, stop_index, distances_sequantial,
                           cluster, cluster, distance_function);
  };

  parallel::parallelize(cluster.size(), fun_parallel, 2);
  parallel::sequential(cluster.size(), fun_sequential, 1);

  EXPECT_TRUE(distances_parallel.isApprox(distances_sequantial));
}

TEST_F(ParallelTest, ReducedEqFullSlice) {
  // auto distance_function = distance::dvstar;

  auto cluster = cluster::get_cluster<container::VLMC_vector>(first_directory, 1);

  matrix_t distances_parallel{cluster.size(), cluster.size()};
  matrix_t distances_sequantial{cluster.size(), cluster.size()};

  auto fun_parallel = [&](size_t start_index, size_t stop_index) {
    calculate::calculate_reduced_slice<container::VLMC_vector>(start_index, stop_index, distances_parallel,
                           cluster, cluster, distance_function);
  };

  auto fun_sequential = [&](size_t start_index, size_t stop_index) {
    calculate::calculate_full_slice<container::VLMC_vector>(start_index, stop_index, distances_sequantial,
                           cluster, cluster, distance_function);
  };

  parallel::parallelize(cluster.size(), fun_parallel, 2);
  parallel::sequential(cluster.size(), fun_sequential, 1);

  size_t y_bound = 1;
  for (size_t i = 1; i < cluster.size(); i++)
  {
    for (size_t j = y_bound; j < y_bound; j++)
    {
      EXPECT_DOUBLE_EQ(distances_parallel(i,j), distances_sequantial(i,j));
    }
    y_bound++;
  }
}