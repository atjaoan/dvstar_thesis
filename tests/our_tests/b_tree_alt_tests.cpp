#pragma once 
#include <gtest/gtest.h>

#include <functional>
#include <filesystem>
#include <cstdlib>
#include "vlmc_from_kmers/kmer.hpp"
#include "b_tree_alt.hpp"

class BTreeTest : public ::testing::Test {
protected:
  void SetUp() override {}
};
TEST_F(BTreeTest, BTreeSearchKmer){
  std::vector<container::RI_Kmer> tmp {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  auto arr = new array::B_Tree(tmp);
  for(int i = 1; i < 16; i++){
    std::cout << arr->a[i].integer_rep << "\n";
  }
}
TEST_F(BTreeTest, SearchKmer){
  std::vector<container::RI_Kmer> tmp {1,2,3};
  auto arr = new array::B_Tree(tmp);
  auto kmer1 = container::RI_Kmer(1);
  auto kmer2 = container::RI_Kmer(2);
  auto kmer3 = container::RI_Kmer(3);
  EXPECT_EQ(kmer1, arr->get_from_array(1));
  EXPECT_EQ(kmer2, arr->get_from_array(2));
  EXPECT_EQ(kmer3, arr->get_from_array(3));
}