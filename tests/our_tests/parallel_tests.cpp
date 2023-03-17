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

  double error_tolerance = 1E-7;
};

TEST_F(ParallelTest, SequentialEqParallel) {
  // auto distance_function = distance::dvstar;
  auto cluster = cluster::get_cluster<container::VLMC_vector>(first_directory, 1, 0);

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

  auto cluster = cluster::get_cluster<container::VLMC_vector>(first_directory, 1, 0);

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

TEST_F(ParallelTest, FullComparisonCheck) {
  auto nr_cores = 1; 
  while(nr_cores <= 8){
    // Sorted Vector Implementation
    auto left_cluster_s = cluster::get_cluster<container::VLMC_sorted_vector>(first_directory, 1, background_order);
    auto right_cluster_s = cluster::get_cluster<container::VLMC_sorted_vector>(first_directory, 1, background_order);
    matrix_t distances_sorted_vector = calculate::calculate_distances<container::VLMC_sorted_vector>(left_cluster_s, right_cluster_s, distance_function, nr_cores);

    // Kmer major implementation
    auto left_cluster_k = cluster::get_kmer_cluster(first_directory, background_order);
    auto right_cluster_k = cluster::get_kmer_cluster(first_directory, background_order);
    matrix_t distances_k_major = calculate::calculate_distance_major(left_cluster_k, right_cluster_k, nr_cores); 

    // Dvstar Original implementation 
    matrix_t distances_org_dvstar{distances_sorted_vector.cols(), distances_sorted_vector.rows()};
    int x = 0;
    for (const auto& dir_entry_x : recursive_directory_iterator(first_directory)) {
      int y = 0; 
      for (const auto& dir_entry_y : recursive_directory_iterator(first_directory)) {
        distances_org_dvstar(x,y) = vlmc::dvstar(dir_entry_x, dir_entry_y, background_order);
        y++;
      }
      x++;
    }
    for (int x = 0; x < distances_sorted_vector.cols(); x++){
      for (int y = 0; y < distances_sorted_vector.rows(); y++){
        if (x==y){
          EXPECT_NEAR(0.0, distances_sorted_vector(x,y), error_tolerance);
          EXPECT_NEAR(0.0, distances_k_major(x,y), error_tolerance);
        } else { 
          EXPECT_NEAR(distances_org_dvstar(x,y), distances_sorted_vector(x,y), error_tolerance);
          EXPECT_NEAR(distances_sorted_vector(x,y), distances_k_major(x,y), error_tolerance); 
        }
      }
    }
    nr_cores = nr_cores * 2; 
  }
}