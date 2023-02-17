#pragma once 
#include <gtest/gtest.h>

#include <functional>
#include <filesystem>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "vlmc_from_kmers/kmer.hpp"
#include "read_in_kmer.hpp"
#include "../read_helper.hpp"

class RIKmerTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

TEST_F(RIKmerTest, EmptyConstructor) {
  container::RI_Kmer kmer{};
  EXPECT_EQ(kmer.integer_rep, 0);
  EXPECT_DOUBLE_EQ(kmer.next_char_prob[0], 0.0);
  EXPECT_DOUBLE_EQ(kmer.next_char_prob[1], 0.0);
  EXPECT_DOUBLE_EQ(kmer.next_char_prob[2], 0.0);
  EXPECT_DOUBLE_EQ(kmer.next_char_prob[3], 0.0);
  EXPECT_TRUE(kmer.is_null);
}
TEST_F(RIKmerTest, KmerConstructorProb) {
  vlmc::VLMCKmer old_kmer{5, 10, {5, 10, 7, 8}};
  double sum_children = 30.0+4.0;
  
  container::RI_Kmer kmer{old_kmer};
  for (size_t i = 0; i < 4; i++){
    EXPECT_DOUBLE_EQ(kmer.next_char_prob[i], (old_kmer.next_symbol_counts[i] + 1)/sum_children);
  }
  EXPECT_FALSE(kmer.is_null);
}

TEST_F(RIKmerTest, KmerConstructorIntRep) {
  std::string kmer_string{"A"};
  auto old_kmer = create_kmer(kmer_string);
  
  container::RI_Kmer kmer{old_kmer};
  EXPECT_EQ(kmer.integer_rep, kmer.get_index_rep(old_kmer));
  EXPECT_FALSE(kmer.is_null);
}

TEST_F(RIKmerTest, KmerConstructorBitRep) {
  std::string kmer_string{"AA"};
  auto old_kmer = create_kmer(kmer_string);
  
  container::RI_Kmer kmer{old_kmer};
  for (size_t i = 0; i < 4; i++){
    EXPECT_EQ(kmer.bit_representation[i], old_kmer.kmer_data[i]);
  }
  EXPECT_FALSE(kmer.is_null);
}
TEST_F(RIKmerTest, KmerBackgroundRep1) {
  std::string kmer_string{"A"};
  auto old_kmer = create_kmer(kmer_string);
  container::RI_Kmer kmer{old_kmer};
  EXPECT_EQ(kmer.background_order_index(kmer.integer_rep, 0), 0);
}
TEST_F(RIKmerTest, KmerBackgroundRep2) {
  std::string kmer_string{"A"};
  auto old_kmer = create_kmer(kmer_string);
  container::RI_Kmer kmer{old_kmer};
  EXPECT_EQ(kmer.background_order_index(kmer.integer_rep, 1), 1);
}
TEST_F(RIKmerTest, KmerBackgroundRep3) {
  std::string kmer_string{"ATT"};
  auto old_kmer = create_kmer(kmer_string);
  container::RI_Kmer kmer{old_kmer};
  EXPECT_EQ(kmer.background_order_index(kmer.integer_rep, 1), 4);
}
TEST_F(RIKmerTest, KmerBackgroundRep4) {
  std::string kmer_string{"ATT"};
  auto old_kmer = create_kmer(kmer_string);
  container::RI_Kmer kmer{old_kmer};
  EXPECT_EQ(kmer.background_order_index(kmer.integer_rep, 2), 20);
}
TEST_F(RIKmerTest, KmerBackgroundRep5) {
  std::string kmer_string{"AAG"};
  auto old_kmer = create_kmer(kmer_string);
  container::RI_Kmer kmer{old_kmer};
  EXPECT_EQ(kmer.background_order_index(kmer.integer_rep, 2), 7);
}