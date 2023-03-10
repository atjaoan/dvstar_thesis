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
  std::filesystem::path path_to_vlmc{"../data/test_VLMCs/sequences_1.bintree"};
};


TEST_F(GetClusterTest, ClusterGetWithVlmcVector) {
  auto container = cluster::get_cluster<container::VLMC_vector>(path_to_bintrees, 1, 0);

  //std::cout << container.get(0).get(0).to_string() << std::endl;
  EXPECT_GT(container.size(), 0);
  EXPECT_GT(container.get(0).size(), 0);
}

TEST_F(GetClusterTest, ClusterGetWithVlmcMultiVector) {
  auto container = cluster::get_cluster<container::VLMC_Indexing>(path_to_bintrees, 1, 0);  
  
  EXPECT_GT(container.size(), 0);
  EXPECT_EQ(container.get(0).get(1).integer_rep, 1);
}

TEST_F(GetClusterTest, ClusterPrettyPrint) {
  container::Kmer_Cluster cluster = cluster::get_kmer_cluster(path_to_bintrees);

  cluster.prettyPrint(); 
}