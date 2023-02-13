#pragma once 
#include <gtest/gtest.h>

#include <functional>
#include <filesystem>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_container.hpp"
#include "cluster_container.hpp"
#include "../read_helper.hpp"

class VlmcContainerTest : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path path_bintree{"../data/test_VLMCs/sequences_1.bintree"};
  //container::VLMC_vector vec_from_path {path_bintree};
};

TEST_F(VlmcContainerTest, IndexRep1) {
  container::VLMC_multi_vector container{};
  std::string kmer_string{"A"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(container.get_index_rep(created), 1);
}

TEST_F(VlmcContainerTest, IndexRep2) {
  container::VLMC_multi_vector container{};
  std::string kmer_string{"T"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(container.get_index_rep(created), 4);
}
TEST_F(VlmcContainerTest, IndexRep3) {
  container::VLMC_multi_vector container{};
  std::string kmer_string{"AA"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(container.get_index_rep(created), 5);
}
TEST_F(VlmcContainerTest, IndexRep4) {
  container::VLMC_multi_vector container{};
  std::string kmer_string{"AT"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(container.get_index_rep(created), 8);
}
TEST_F(VlmcContainerTest, IndexRep5) {
  container::VLMC_multi_vector container{};
  std::string kmer_string{"CC"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(container.get_index_rep(created), 10);
}
TEST_F(VlmcContainerTest, IndexRep6) {
  container::VLMC_multi_vector container{};
  std::string kmer_string{"TG"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(container.get_index_rep(created), 19);
}
TEST_F(VlmcContainerTest, IndexRep7) {
  container::VLMC_multi_vector container{};
  std::string kmer_string{"AAA"};
  auto created = create_kmer(kmer_string);
  EXPECT_EQ(container.get_index_rep(created), 21);
}

TEST_F(VlmcContainerTest, AddToMultiVector) {
  std::ifstream ifs(path_bintree, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);
  container::VLMC_multi_vector container{};
  
  Kmer kmer{};

  if (ifs.peek() != EOF){
    archive(kmer);
    container.push(kmer);
  }
  ifs.close();
  
  int kmer_index = container.get_index_rep(kmer);
  std::cout << kmer_index << " : " << kmer.to_string() << std::endl;
  EXPECT_EQ(container.get(kmer_index).to_string(), kmer.to_string());
}

TEST_F(VlmcContainerTest, IndexAdd1) {
  container::VLMC_multi_vector container{};
  std::string kmer_string{"A"};
  std::cout << "1" << std::endl;
  Kmer created = create_kmer(kmer_string);
  std::cout << "2" << std::endl;
  //std::cout << container.size() << std::endl;
  std::cout << "3" << std::endl;
  container.push(created);
  std::cout << "end" << std::endl;
  //EXPECT_EQ(container.get(1).to_string(), "A");
}
