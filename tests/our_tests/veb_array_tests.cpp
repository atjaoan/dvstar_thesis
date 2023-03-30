#pragma once 
#include <gtest/gtest.h>

#include <functional>
#include <filesystem>
#include <cstdlib>
#include "vlmc_from_kmers/kmer.hpp"
#include "veb_array.hpp"

class VebArrayTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

TEST_F(VebArrayTest, Indicies) {
  veb::Veb_array tree(16);
  int k = 4; // height
	int pow_k = tree.power_of_two(k); // = 4^2 = 16

  EXPECT_EQ(5, tree.get_index(1, k));
  EXPECT_EQ(4, tree.get_index(2, k));
  EXPECT_EQ(6, tree.get_index(3, k));
  EXPECT_EQ(2, tree.get_index(4, k));
  EXPECT_EQ(8, tree.get_index(5, k));
  EXPECT_EQ(7, tree.get_index(6, k));
  EXPECT_EQ(9, tree.get_index(7, k));
  EXPECT_EQ(1, tree.get_index(8, k));
  EXPECT_EQ(11, tree.get_index(9, k));
  EXPECT_EQ(10, tree.get_index(10, k));
  EXPECT_EQ(12, tree.get_index(11, k));
  EXPECT_EQ(3, tree.get_index(12, k));
  EXPECT_EQ(14, tree.get_index(13, k));
  EXPECT_EQ(13, tree.get_index(14, k));
  EXPECT_EQ(15, tree.get_index(15, k));
}

TEST_F(VebArrayTest, Size10TreeInsert){
  veb::Veb_array tree(10);
  int power = tree.power_of_two((int)(ceil(log((float)10)/log(2.0))));

  int height = tree.height;

  int idx = tree.get_index(5, tree.height);
  tree.container[idx] = 5;
  for(int i = 0; i < tree.size; i++){
    if(i != 8)
      EXPECT_EQ(0, tree.container[i]);
  }
  EXPECT_EQ(5, tree.container[8]);
  EXPECT_EQ(16, power);
  EXPECT_EQ(4, tree.height);
}

TEST_F(VebArrayTest, InsertFn){
  veb::Veb_array tree(10);
  int power = tree.power_of_two((int)(ceil(log((float)10)/log(2.0))));

  int height = tree.height;

  tree.insert(5);
  for(int i = 0; i < tree.size; i++){
    if(i != 8)
      EXPECT_EQ(0, tree.container[i]);
  }
  EXPECT_EQ(5, tree.container[8]);
  EXPECT_EQ(16, power);
  EXPECT_EQ(4, tree.height);
}

TEST_F(VebArrayTest, VebSearch){
  veb::Veb_array tree(16);
  tree.insert(5);
  tree.insert(1);
  tree.insert(2);
  EXPECT_EQ(5, tree.get_elem(5));
  EXPECT_EQ(1, tree.get_elem(1));
  EXPECT_EQ(2, tree.get_elem(2));
  //tree.insert(9);
  //std::cout << tree.veb_search(tree.container.data(), tree.height, 1) << "\n";
  //std::cout << tree.veb_search(tree.container.data(), tree.size, 5) << "\n";
  //std::cout << tree.veb_search(tree.container.data(), tree.size, 2) << "\n";
  //EXPECT_EQ(5, tree.veb_search(tree.container.data(), tree.size, 5));
  //EXPECT_EQ(1, tree.veb_search(tree.container.data(), tree.height, 1));
  //EXPECT_EQ(2, tree.veb_search(tree.container.data(), tree.height, 2));
  //EXPECT_EQ(9, tree.veb_search(tree.container.data(), tree.height, 9));
}