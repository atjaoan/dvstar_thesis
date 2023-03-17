#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

#include "vlmc_container.hpp"
#include "cluster_container.hpp"

using vlmc_c = container::VLMC_vector;
using cluster_c = container::Cluster_Container<vlmc_c>;

class ClusterContainerTest : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_bintree{"../data/test_VLMCs/sequences_1.bintree"};
  std::filesystem::path second_bintree{"../data/test_VLMCs/sequences_2.bintree"};
  std::filesystem::path third_bintree{"../data/test_VLMCs/sequences_3.bintree"};

  vlmc_c first_vlmc{first_bintree};
  vlmc_c second_vlmc{second_bintree};
  vlmc_c third_vlmc{third_bintree};
};


TEST_F(ClusterContainerTest, AddToContainer) {
  cluster_c container {};
  container.push(first_vlmc);
  container.push(second_vlmc);
  container.push(third_vlmc);
  EXPECT_EQ(container.size(), 3);
}

TEST_F(ClusterContainerTest, VlmcSizeNonZeroAfterAddToContainer) {
  cluster_c container {};
  container.push(first_vlmc);
  EXPECT_GT(container.get(0).size(), 0);
}

TEST_F(ClusterContainerTest, KmerContainerGet){
  container::Kmer_Cluster container{};
  container::Kmer_Pair kmer0 = container::Kmer_Pair(container::RI_Kmer(0), 0);
  container.push(kmer0);
  std::vector<container::Kmer_Pair> kmer_bucket = container.get(container.get_bucket(kmer0));
  EXPECT_EQ(kmer0.id, kmer_bucket[0].id);
}