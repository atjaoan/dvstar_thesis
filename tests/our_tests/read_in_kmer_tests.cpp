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

void createStrComb(std::vector<std::string> &strings, std::string str, int size){
  if (size <= 0){
    return; 
  }
  std::string str_A = str; 
  std::string str_C = str; 
  std::string str_G = str;
  std::string str_T = str; 

  str_A.append("A");
  strings.push_back(str_A);
  str_C.append("C");
  strings.push_back(str_C);
  str_G.append("G");
  strings.push_back(str_G);
  str_T.append("T");
  strings.push_back(str_T);

  createStrComb(strings, str_A, size - 1);
  createStrComb(strings, str_C, size - 1);
  createStrComb(strings, str_G, size - 1);
  createStrComb(strings, str_T, size - 1);
}

bool compare (std::pair <int, std::string> &a, std::pair <int, std::string> &b ){
   return a.first < b.first;
}

std::string get_background_context(const std::string &state,
                                   const size_t background_order) {
  if (state.size() <= background_order) {
    return state; // <- This will never happen 
  } else {
    size_t background = state.size() - background_order;
    return state.substr(background);
  }
}

TEST_F(RIKmerTest, showIntRep) {
  container::RI_Kmer ri_kmer{};
  std::vector<std::string> strings{""};
  std::string current_string{""};  
  createStrComb(strings, current_string, 3);

  std::vector<std::pair<int, std::string>> cmb{};

  for (auto str : strings){
    std::string kmer_string{str};
    auto old_kmer = create_kmer(kmer_string);
    container::RI_Kmer kmer{old_kmer};
    cmb.push_back({kmer.integer_rep, str}); 
  }

  sort(cmb.begin(), cmb.end(), compare);

  // for (auto item : cmb) {
  //   std::cout << "Index " << item.first << " -> " << item.second << std::endl; 
  // }

  for (int background_order = 2; background_order < 3; background_order++){
    for (auto item : cmb) {
      auto idx = item.first;
      auto context = get_background_context(item.second, background_order);
      std::string kmer_string{context};
      auto old_kmer = create_kmer(kmer_string);
      container::RI_Kmer kmer{old_kmer};
      auto context_idx = kmer.integer_rep;
      auto our_context_idx = ri_kmer.background_order_index(idx, background_order);
      std::cout << "Index " << idx << " " << item.second << std::endl; 
      EXPECT_EQ(context_idx, our_context_idx);
    }
  }

}