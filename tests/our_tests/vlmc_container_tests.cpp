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
#include "read_in_kmer.hpp"

class VlmcContainerTest : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path path_bintree{"../data/test_VLMCs/sequences_1.bintree"};
  //container::VLMC_vector vec_from_path {path_bintree};
};

TEST_F(VlmcContainerTest, IndexRep1) {
  container::Index_by_value container{};
  std::string kmer_string{"A"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 1);
}

TEST_F(VlmcContainerTest, IndexRep2) {
  container::Index_by_value container{};
  std::string kmer_string{"T"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 4);
}
TEST_F(VlmcContainerTest, IndexRep3) {
  container::Index_by_value container{};
  std::string kmer_string{"AA"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 5);
}
TEST_F(VlmcContainerTest, IndexRep4) {
  container::Index_by_value container{};
  std::string kmer_string{"AT"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 8);
}
TEST_F(VlmcContainerTest, IndexRep5) {
  container::Index_by_value container{};
  std::string kmer_string{"CC"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 10);
}
TEST_F(VlmcContainerTest, IndexRep6) {
  container::Index_by_value container{};
  std::string kmer_string{"TG"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 19);
}
TEST_F(VlmcContainerTest, IndexRep7) {
  container::Index_by_value container{};
  std::string kmer_string{"AAA"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 21);
}

TEST_F(VlmcContainerTest, AddReadInKmerToIndexByValue) {
  container::Index_by_value container{};
  std::string kmer_string{"A"};
  Kmer created = create_kmer(kmer_string);
  container::RI_Kmer in_kmer{created};
  container.push(in_kmer);
  EXPECT_EQ(container.get(1), in_kmer);
}

TEST_F(VlmcContainerTest, AddManyReadInKmerToIndexByValue) {
  container::Index_by_value container{};
  std::vector<std::string> kmer_strings {"A", "C", "G", "T", "AC"};
  std::vector<container::RI_Kmer> in_kmers {};
  for(auto s : kmer_strings){
    Kmer created = create_kmer(s);
    container::RI_Kmer new_kmer{created};
    in_kmers.push_back(new_kmer);
    container.push(new_kmer);
  } 

  for(auto kmer : in_kmers){
    EXPECT_EQ(kmer, container.get(kmer.integer_rep));
  }
}
/*
TEST_F(VlmcContainerTest, UseDirectoryIndexByValue) {
  container::Index_by_value container{path_bintree};
  container::VLMC_multi_vector multi_vec{path_bintree};
  for(size_t i = 1; i <= multi_vec.get_max_kmer_index(); i++){
    if(multi_vec.get_index_rep(multi_vec.get(i)) == 0) continue;
    EXPECT_EQ(container.get(i).integer_rep, i);
  }
}
*/
TEST_F(VlmcContainerTest, TimingTest) {
  int items = 6561;
  container::VLMC_sorted_vector tree{};
  auto null_kmer = container::RI_Kmer(-1);
  std::chrono::steady_clock::time_point begin_insert = std::chrono::steady_clock::now();
  for (size_t i = 0; i < items; i++){
    auto kmer = container::RI_Kmer{i};
    tree.push(kmer);
  }
  std::chrono::steady_clock::time_point end_insert = std::chrono::steady_clock::now();
  auto insert_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_insert - begin_insert).count();
  std::cout << "Time insert : " << insert_time << " [nano sec]" << std::endl;
  std::cout << "sec / items : " << insert_time / items << " [nano sec]" << std::endl;

  auto begin_find = std::chrono::steady_clock::now();
  for (size_t i = 0; i < items; i++){
    tree.find(i);
  }
  auto end_find = std::chrono::steady_clock::now();
  auto find_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_find - begin_find).count();
  std::cout << "Time find : " << find_time << " [nano sec]" << std::endl;
  std::cout << "sec / items : " << find_time / items << " [nano sec]" << std::endl;
}
/*
TEST_F(VlmcContainerTest, VLMC_IndexingCaching) {
  container::Index_by_value container{};
  container::RI_Kmer in_kmer{created};
  container.push(in_kmer);
  EXPECT_EQ(container.get(1), in_kmer);
}
*/