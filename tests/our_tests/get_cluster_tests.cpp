#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

//Our files 
#include "vlmc_container.hpp"
#include "get_cluster.hpp"

using VLMC_vector = container::VLMC_vector;
using Index_by_value = container::Index_by_value;

class GetClusterTest : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path path_to_bintrees{"../data/test_VLMCs"};
  std::filesystem::path path_to_vlmc{"../data/test_VLMCs/sequences_1.bintree"};
};


TEST_F(GetClusterTest, ClusterGetWithVlmcVector) {
  container::Cluster_vector container{};
  cluster::get_cluster<VLMC_vector>(path_to_bintrees, container);

  //std::cout << container.get(0).get(0).to_string() << std::endl;
  EXPECT_GT(container.size(), 0);
  EXPECT_GT(container.get(0).size(), 0);
}

TEST_F(GetClusterTest, ClusterGetWithVlmcMultiVector) {
  container::Cluster_vector container{};
  cluster::get_cluster<Index_by_value>(path_to_bintrees, container); 

  EXPECT_GT(container.size(), 0);
  EXPECT_EQ(container.get(0).get(1).integer_rep, 1);
}