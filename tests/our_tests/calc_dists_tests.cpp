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

//Original implementation files
#include <kmc_file.h>
#include "vlmc_from_kmers/dvstar.hpp"
#include "vlmc_from_kmers/build_vlmc.hpp"

using matrix_t  = Eigen::MatrixXd;
using vlmc_c = container::VLMC_vector;
using vlmc_t = container::VLMC_Container;
using cluster_c = container::Cluster_vector;

class CalcDistsTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_bintree{"../data/test_VLMCs/sequences_1.bintree"};
  std::filesystem::path second_bintree{"../data/test_VLMCs/sequences_2.bintree"};
  std::filesystem::path third_bintree{"../data/test_VLMCs/sequences_3.bintree"};

  vlmc_c first_vlmc{first_bintree};
  vlmc_c second_vlmc{second_bintree};
  vlmc_c third_vlmc{third_bintree};

  std::function<double(vlmc_t &, vlmc_t &)> dist_func = distance::dvstar;
};


TEST_F(CalcDistsTests, Identity) {
  cluster_c cluster_left{};
  cluster_c cluster_right{};
  cluster_left.push(std::make_shared<vlmc_c>(first_vlmc));
  cluster_left.push(std::make_shared<vlmc_c>(second_vlmc));
  cluster_right.push(std::make_shared<vlmc_c>(first_vlmc));
  cluster_right.push(std::make_shared<vlmc_c>(second_vlmc));

  matrix_t distances = calculate::calculate_distances(cluster_left, dist_func, 1);
  EXPECT_DOUBLE_EQ(distances(0,0), 0.0);
  EXPECT_DOUBLE_EQ(distances(1,1), 0.0);
}
