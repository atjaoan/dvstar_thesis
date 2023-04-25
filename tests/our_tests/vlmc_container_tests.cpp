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
  std::filesystem::path path_corner_case{"../data/medium_test/sequences_951.fasta.bintree"};
  //container::VLMC_vector vec_from_path {path_bintree};
};

TEST_F(VlmcContainerTest, IndexRep1) {
  container::VLMC_sorted_vector container{};
  std::string kmer_string{"A"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 1);
}

TEST_F(VlmcContainerTest, IndexRep2) {
  container::VLMC_sorted_vector container{};
  std::string kmer_string{"T"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 4);
}
TEST_F(VlmcContainerTest, IndexRep3) {
  container::VLMC_sorted_vector container{};
  std::string kmer_string{"AA"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 5);
}
TEST_F(VlmcContainerTest, IndexRep4) {
  container::VLMC_sorted_vector container{};
  std::string kmer_string{"AT"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 8);
}
TEST_F(VlmcContainerTest, IndexRep5) {
  container::VLMC_sorted_vector container{};
  std::string kmer_string{"CC"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 10);
}
TEST_F(VlmcContainerTest, IndexRep6) {
  container::VLMC_sorted_vector container{};
  std::string kmer_string{"TG"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 19);
}
TEST_F(VlmcContainerTest, IndexRep7) {
  container::VLMC_sorted_vector container{};
  std::string kmer_string{"AAA"};
  auto created = create_kmer(kmer_string);
  container::RI_Kmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 21);
}