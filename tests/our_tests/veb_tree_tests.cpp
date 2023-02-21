#pragma once 
#include <gtest/gtest.h>

#include <functional>
#include <filesystem>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "vlmc_from_kmers/kmer.hpp"
#include "veb_tree.hpp"
#include "vlmc_container.hpp"
#include "cluster_container.hpp"
#include "../read_helper.hpp"

class VebTreeTest : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path path_bintree{"../data/test_VLMCs/sequences_1.bintree"};
  container::RI_Kmer max_kmer = container::RI_Kmer(INT_MAX);
  container::RI_Kmer min_kmer = container::RI_Kmer(INT_MIN);
};
TEST_F(VebTreeTest, SizedConstructor) {
  veb::Veb_tree tree{16};
  EXPECT_EQ(tree.min, max_kmer);
  EXPECT_EQ(tree.max, min_kmer);
  EXPECT_EQ(tree.size, 16);
  EXPECT_NE(nullptr, tree.summary);
  EXPECT_EQ(nullptr, tree.trees[0]);
  EXPECT_TRUE(tree.is_empty);
}
TEST_F(VebTreeTest, InsertOne) {
  veb::Veb_tree tree{16};
  auto kmer =  container::RI_Kmer{1};
  veb::insert(tree, kmer);
  EXPECT_EQ(tree.min, kmer);
  EXPECT_EQ(tree.max, kmer);
  EXPECT_EQ(tree.size, 16);
  EXPECT_EQ(nullptr, tree.trees[0]);
  EXPECT_NE(nullptr, tree.summary);
}
TEST_F(VebTreeTest, SummaryTreeAfterTwoInsert) {
  veb::Veb_tree tree{16};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer2 = container::RI_Kmer{4};
  //std::cout << "Kmer two has integer rep : " << kmer2.integer_rep << std::endl;
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer2);
  //std::cout << "Kmer one has integer rep : " << kmer1.integer_rep << std::endl;
  //std::cout << "Kmer two has integer rep : " << kmer2.integer_rep << std::endl;
  EXPECT_NE(nullptr, tree.summary);
  EXPECT_EQ(tree.summary->size, 4);
  EXPECT_EQ(tree.max, kmer2);
  EXPECT_EQ(tree.min, kmer1);
}

TEST_F(VebTreeTest, SubTreeAfterTwoInsert) {
  veb::Veb_tree tree{16};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer0 = container::RI_Kmer{0};
  veb::insert(tree, kmer1);
  EXPECT_EQ(kmer1, tree.min);
  EXPECT_EQ(kmer1, tree.max);
  veb::insert(tree, kmer0);
  EXPECT_EQ(kmer0, tree.min);
  EXPECT_NE(nullptr, tree.summary);
  EXPECT_EQ(kmer0, tree.summary->min);
  EXPECT_EQ(kmer0, tree.summary->max);
  EXPECT_EQ(tree.trees[0]->min, kmer1);
  EXPECT_EQ(tree.trees[0]->max, kmer1);
  EXPECT_EQ(tree.min, kmer0);
}
TEST_F(VebTreeTest, FindOnOneInsert) {
  veb::Veb_tree tree{16};
  auto kmer1 = container::RI_Kmer{1};
  veb::insert(tree, kmer1);
  EXPECT_EQ(kmer1, veb::find(tree, kmer1));
}

TEST_F(VebTreeTest, FindOnTwoInsert) {
  veb::Veb_tree tree{16};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer0 = container::RI_Kmer{0};
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer0);
  EXPECT_EQ(kmer1, veb::find(tree, kmer1));
  EXPECT_EQ(kmer0, veb::find(tree, kmer0));
}
/*
TEST_F(VebTreeTest, Insert16Find16) {
  std::vector<container::RI_Kmer> kmers {};
  for (size_t i = 0; i < 16; i++){
    auto kmer1 = container::RI_Kmer{i};
    kmers.push_back(container::RI_Kmer{i});
  }
  
  veb::Veb_tree tree{16};
  for (auto kmer : kmers){
    veb::insert(tree, kmer);
  }

  for (int i = 0; i < 16; i++){
    EXPECT_EQ(kmers[i], veb::find(tree, kmers[i]));
  }
}
*/
TEST_F(VebTreeTest, FindOneSucc) {
  veb::Veb_tree tree{16};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer2 = container::RI_Kmer{2};
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer2);
  //std::cout << "max index is " << tree.max.integer_rep << std::endl;
  EXPECT_EQ(kmer2, veb::succ(tree, kmer1));
}

TEST_F(VebTreeTest, Find4) {
  veb::Veb_tree tree{16};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer2 = container::RI_Kmer{2};
  auto kmer3 = container::RI_Kmer{3};
  auto kmer4 = container::RI_Kmer{4};
  auto kmer5 = container::RI_Kmer{5};
  std::cout << "1" << std::endl;
  veb::insert(tree, kmer1);
  std::cout << "2" << std::endl;
  veb::insert(tree, kmer2);
  std::cout << "3" << std::endl;
  veb::insert(tree, kmer3);
  std::cout << "4" << std::endl;
  veb::insert(tree, kmer4);
  std::cout << "5" << std::endl;
  veb::insert(tree, kmer5);
  std::cout << "succs" << std::endl;
  std::cout << "1" << std::endl;
  EXPECT_EQ(kmer2, veb::succ(tree, kmer1));
  std::cout << "2" << std::endl;
  EXPECT_EQ(kmer3, veb::succ(tree, kmer2));
  std::cout << "3" << std::endl;
  EXPECT_EQ(kmer4, veb::succ(tree, kmer3));
  std::cout << "4" << std::endl;
  EXPECT_EQ(kmer5, veb::succ(tree, kmer4));
}

/*
*/