#pragma once 

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

//Our files 
#include "vlmc_container.hpp"
#include "cluster_container.hpp"
#include "calc_dists.hpp"
#include "distances/dvstar.hpp"
#include "get_cluster.hpp"
#include "global_aliases.hpp"

//Original implementation files
#include <kmc_file.h>
#include "vlmc_from_kmers/dvstar.hpp"
#include "vlmc_from_kmers/build_vlmc.hpp"

using cluster_c = container::Cluster_Container<container::VLMC_vector>;
extern const out_t error_tolerance;

class CalcDistsTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_bintree{"../data/test_VLMCs/sequences_1.bintree"};
  std::filesystem::path second_bintree{"../data/test_VLMCs/sequences_2.bintree"};
  std::filesystem::path third_bintree{"../data/test_VLMCs/sequences_3.bintree"};

  std::filesystem::path path_to_bintrees{"../data/test_VLMCs"};

  container::VLMC_vector first_vlmc{first_bintree};
  container::VLMC_vector second_vlmc{second_bintree};
  container::VLMC_vector third_vlmc{third_bintree};

  size_t background_order = 0;

  std::function<out_t(container::VLMC_vector &, container::VLMC_vector &)> dist_func = [&](auto &left, auto &right) {
      return distance::dvstar<container::VLMC_vector>(left, right, background_order);
  };

  std::function<out_t(container::VLMC_sorted_vector &, container::VLMC_sorted_vector &)> dist_func_sv = [&](auto &left, auto &right) {
      return distance::dvstar<container::VLMC_sorted_vector>(left, right, background_order);
  };

  std::function<out_t(container::VLMC_Indexing &, container::VLMC_Indexing &)> dist_func_mv = [&](auto &left, auto &right) {
      return distance::dvstar<container::VLMC_Indexing>(left, right, background_order);
  };

  std::function<out_t(container::VLMC_B_tree &, container::VLMC_B_tree &)> dist_func_b = [&](auto &left, auto &right) {
      return distance::dvstar<container::VLMC_B_tree>(left, right, background_order);
  };

  std::function<out_t(container::VLMC_hashmap &, container::VLMC_hashmap &)> dist_func_h = [&](auto &left, auto &right) {
      return distance::dvstar<container::VLMC_hashmap>(left, right, background_order);
  };

  std::function<out_t(container::VLMC_Combo &, container::VLMC_Combo &)> dist_func_c = [&](auto &left, auto &right) {
      return distance::dvstar<container::VLMC_Combo>(left, right, background_order);
  };

  std::function<out_t(container::VLMC_Veb &, container::VLMC_Veb &)> dist_func_veb = [&](auto &left, auto &right) {
      return distance::dvstar<container::VLMC_Veb>(left, right, background_order);
  };

  std::function<out_t(container::VLMC_Set &, container::VLMC_Set &)> dist_func_set = [&](auto &left, auto &right) {
      return distance::dvstar<container::VLMC_Set>(left, right, background_order);
  };
};

cluster_c create_cluster(container::VLMC_vector &vlmc, size_t size) {
  cluster_c cluster{};
  for (size_t i = 0; i < size; i++){
    cluster.push(vlmc);
  }
  return cluster; 
}

// Tests when comparing two directories 
TEST_F(CalcDistsTests, SizeTwoDir) {
  for (size_t x = 1; x < 5; x++){
    for (size_t y = 1; y < 5; y++){
      cluster_c cluster_left = create_cluster(first_vlmc, x);
      cluster_c cluster_right = create_cluster(second_vlmc, y);

      matrix_t distances = calculate::calculate_distances<container::VLMC_vector>(cluster_left, cluster_right, dist_func, 1);
      EXPECT_EQ(distances.size(), x*y);
      EXPECT_EQ(distances.rows(), x);
      EXPECT_EQ(distances.cols(), y);
    }
  }
}

TEST_F(CalcDistsTests, AllValsTwoDir) {
  for (size_t x = 1; x < 2; x++){
    cluster_c cluster_left = create_cluster(first_vlmc, x);
    cluster_c cluster_right = create_cluster(second_vlmc, x);

    matrix_t distances = calculate::calculate_distances<container::VLMC_vector>(cluster_left, cluster_right, dist_func, 1);

    for (size_t row = 0; row < distances.rows(); row++){
      for (size_t col = 0; col < distances.cols(); col++){
        EXPECT_GT(distances(row,col), 0);
      }
    }
  }
}

TEST_F(CalcDistsTests, ValueCheckTwoDir){
  // Vector Implementation
  auto left_cluster_v = cluster::get_cluster<container::VLMC_vector>(path_to_bintrees, 1, background_order);
  auto right_cluster_v = cluster::get_cluster<container::VLMC_vector>(path_to_bintrees, 1, background_order);
  matrix_t distances_vector = calculate::calculate_distances<container::VLMC_vector>(left_cluster_v, right_cluster_v, dist_func, 1);

  // Multi Vector Implementation
  auto left_cluster_mv = cluster::get_cluster<container::VLMC_Indexing>(path_to_bintrees, 1, background_order);
  auto right_cluster_mv = cluster::get_cluster<container::VLMC_Indexing>(path_to_bintrees, 1, background_order);
  matrix_t distances_indexing = calculate::calculate_distances<container::VLMC_Indexing>(left_cluster_mv, right_cluster_mv, dist_func_mv, 1);
 
  // Sorted Vector Implementation
  auto left_cluster_s = cluster::get_cluster<container::VLMC_sorted_vector>(path_to_bintrees, 1, background_order);
  auto right_cluster_s = cluster::get_cluster<container::VLMC_sorted_vector>(path_to_bintrees, 1, background_order);
  matrix_t distances_sorted_vector = calculate::calculate_distances<container::VLMC_sorted_vector>(left_cluster_s, right_cluster_s, dist_func_sv, 1); 

  // B-tree Implementation
  auto left_cluster_b = cluster::get_cluster<container::VLMC_B_tree>(path_to_bintrees, 1, background_order);
  auto right_cluster_b = cluster::get_cluster<container::VLMC_B_tree>(path_to_bintrees, 1, background_order);
  matrix_t distances_b_tree = calculate::calculate_distances<container::VLMC_B_tree>(left_cluster_b, right_cluster_b, dist_func_b, 1); 

  // HashMap Implementation
  auto left_cluster_h = cluster::get_cluster<container::VLMC_hashmap>(path_to_bintrees, 1, background_order);
  auto right_cluster_h = cluster::get_cluster<container::VLMC_hashmap>(path_to_bintrees, 1, background_order);
  matrix_t distances_hashmap = calculate::calculate_distances<container::VLMC_hashmap>(left_cluster_h, right_cluster_h, dist_func_h, 1); 

  // Combination of Index by value & sorted vector Implementation
  auto left_cluster_c = cluster::get_cluster<container::VLMC_Combo>(path_to_bintrees, 1, background_order);
  auto right_cluster_c = cluster::get_cluster<container::VLMC_Combo>(path_to_bintrees, 1, background_order);
  matrix_t distances_combo = calculate::calculate_distances<container::VLMC_Combo>(left_cluster_c, right_cluster_c, dist_func_c, 1); 

  // Kmer major implementation
  auto left_cluster_k = cluster::get_kmer_cluster(path_to_bintrees, background_order);
  auto right_cluster_k = cluster::get_kmer_cluster(path_to_bintrees, background_order);
  matrix_t distances_k_major = calculate::calculate_distance_major(left_cluster_k, right_cluster_k, 1);  

  // Veb Implementation - Currently too slow
  auto left_cluster_veb = cluster::get_cluster<container::VLMC_Veb>(path_to_bintrees, 1, background_order);
  auto right_cluster_veb = cluster::get_cluster<container::VLMC_Veb>(path_to_bintrees, 1, background_order);
  matrix_t distances_veb = calculate::calculate_distances<container::VLMC_Veb>(left_cluster_veb, right_cluster_veb, dist_func_veb, 1); 

  // std::Set
  auto left_cluster_set = cluster::get_cluster<container::VLMC_Set>(path_to_bintrees, 1, background_order);
  auto right_cluster_set = cluster::get_cluster<container::VLMC_Set>(path_to_bintrees, 1, background_order);
  matrix_t distances_set = calculate::calculate_distances<container::VLMC_Set>(left_cluster_set, right_cluster_set, dist_func_set, 1); 

  // Dvstar Original implementation 
  matrix_t distances_org_dvstar{distances_vector.cols(), distances_vector.rows()};
  int x = 0;
  for (const auto& dir_entry_x : recursive_directory_iterator(path_to_bintrees)) {
    int y = 0; 
    for (const auto& dir_entry_y : recursive_directory_iterator(path_to_bintrees)) {
      distances_org_dvstar(x,y) = vlmc::dvstar(dir_entry_x, dir_entry_y, background_order);
      y++;
    }
    x++;
  }
  for (int x = 0; x < distances_vector.cols(); x++){
    for (int y = 0; y < distances_vector.rows(); y++){
      if (x==y){
        // EXPECT_NEAR(0.0, distances_vector(x,y), error_tolerance);
        // EXPECT_NEAR(0.0, distances_indexing(x,y), error_tolerance);
      } else { 
        EXPECT_NEAR(distances_org_dvstar(x,y), distances_vector(x,y), error_tolerance);
        EXPECT_NEAR(distances_vector(x,y), distances_sorted_vector(x,y), error_tolerance);
        EXPECT_NEAR(distances_vector(x,y), distances_indexing(x,y), error_tolerance);
        EXPECT_NEAR(distances_vector(x,y), distances_b_tree(x,y), error_tolerance);
        EXPECT_NEAR(distances_vector(x,y), distances_hashmap(x,y), error_tolerance);
        EXPECT_NEAR(distances_vector(x,y), distances_combo(x,y), error_tolerance);
        // EXPECT_NEAR(distances_vector(x,y), distances_veb(x,y), error_tolerance);
        EXPECT_NEAR(distances_vector(x,y), distances_set(x,y), error_tolerance);
        EXPECT_NEAR(distances_vector(x,y), distances_k_major(x,y), error_tolerance); 
      }
    }
  }
}

TEST_F(CalcDistsTests, ValueCheckOneDir){
  // Vector Implementation
  auto left_cluster_v = cluster::get_cluster<container::VLMC_vector>(path_to_bintrees, 1, background_order);
  auto right_cluster_v = cluster::get_cluster<container::VLMC_vector>(path_to_bintrees, 1, background_order);
  matrix_t distances_vector_two_dirs = calculate::calculate_distances<container::VLMC_vector>(left_cluster_v, right_cluster_v, dist_func, 1);

  auto cluster_v = cluster::get_cluster<container::VLMC_vector>(path_to_bintrees, 1, background_order);

  for (int nr_cores = 1; nr_cores < 9; nr_cores++){
    matrix_t distances_vector = calculate::calculate_distances<container::VLMC_vector>(cluster_v, dist_func, nr_cores);
    for (int x = 0; x < distances_vector.cols(); x++){
      for (int y = 0; y < distances_vector.rows(); y++){
        if (x <= y) {
          EXPECT_NEAR(distances_vector_two_dirs(x,y), distances_vector(x,y), error_tolerance);
        } else {
          EXPECT_NEAR(0.0, distances_vector(x,y), error_tolerance);
        }
      }
    }
  }

}

// Tests for inter-directory comparisons
TEST_F(CalcDistsTests, SizeOneDir) {
  for (size_t x = 1; x < 5; x++){
    cluster_c cluster_left = create_cluster(first_vlmc, x);
    matrix_t distances = calculate::calculate_distances(cluster_left, dist_func, 1);
    EXPECT_EQ(distances.size(), x * x); 
    EXPECT_EQ(distances.rows(), x);
    EXPECT_EQ(distances.cols(), x);
  }
}

TEST_F(CalcDistsTests, AllValsOneDir) {
  for (size_t x = 1; x < 2; x++){
    cluster_c cluster = create_cluster(first_vlmc, x);

    matrix_t distances = calculate::calculate_distances(cluster, dist_func, 1);

    for (size_t row = 0; row < distances.rows(); row++){
      for (size_t col = 1 + row; col < distances.cols(); col++){
        EXPECT_GT(distances(row,col), 0);
      }
    }
  }
}

