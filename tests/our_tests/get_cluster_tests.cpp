#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

//Our files 
#include "vlmc_container.hpp"
#include "get_cluster.hpp"

class GetClusterTest : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path path_to_bintrees{"../data/test_VLMCs"};
};


TEST_F(GetClusterTest, ClusterGetWithVlmcVector) {
  container::Cluster_vector container{};
  cluster::get_cluster<container::VLMC_vector>(path_to_bintrees, container);


  EXPECT_GT(container.size(), 0);
  EXPECT_GT(container.get(0).size(), 0);
}