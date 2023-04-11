#pragma once 
#include <gtest/gtest.h>

#include <functional>
#include <filesystem>
#include <cstdlib>
#include "vlmc_from_kmers/kmer.hpp"
#include "eytzinger_array.hpp"

class EytzArrayTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

TEST_F(EytzArrayTest, SearchKmer){
  std::vector<container::RI_Kmer> tmp {1,2,3};
  auto arr = new array::Ey_array(tmp);
  auto kmer1 = container::RI_Kmer(1);
  auto kmer2 = container::RI_Kmer(2);
  auto kmer3 = container::RI_Kmer(3);
  EXPECT_EQ(kmer1, arr->get_from_array(1));
  EXPECT_EQ(kmer2, arr->get_from_array(2));
  EXPECT_EQ(kmer3, arr->get_from_array(3));
}
TEST_F(EytzArrayTest, Construct){
  std::vector<container::RI_Kmer> tmp {1,2,3, 4, 5, 6, 7,8, 9, 10, 11, 12};
  auto arr = new array::Ey_array(tmp);
  for(int i = 1; i < 16; i++){
    //std::cout << arr->b[i].integer_rep << "\n";
  }
}