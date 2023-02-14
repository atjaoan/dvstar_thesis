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
  container::Cluster_vector container{};
  cluster::get_cluster<container::VLMC_vector>(path_to_bintrees, container);

  //std::cout << container.get(0).get(0).to_string() << std::endl;
  EXPECT_GT(container.size(), 0);
  EXPECT_GT(container.get(0).size(), 0);
}

TEST_F(GetClusterTest, ClusterGetWithVlmcMultiVector) {
  container::Cluster_vector container{};
  cluster::get_cluster<container::VLMC_multi_vector>(path_to_bintrees, container); 

  std::ifstream ifs(path_to_vlmc, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);
  container::VLMC_multi_vector vlmc_multi_vec{};
  
  Kmer kmer{};

  if (ifs.peek() != EOF){
    archive(kmer);
  }
  ifs.close();

  int index = vlmc_multi_vec.get_index_rep(kmer);
  container::VLMC_Container vlmc = container.get(0);
  EXPECT_GT(container.size(), 0);
  EXPECT_EQ(vlmc.get(index).to_string(), kmer.to_string());
}