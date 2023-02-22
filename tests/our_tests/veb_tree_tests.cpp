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
  container::RI_Kmer max_kmer = container::RI_Kmer(-1);
  container::RI_Kmer min_kmer = container::RI_Kmer(INT_MAX);
};

TEST_F(VebTreeTest, SizedConstructor) {
  veb::Veb_tree tree{16};
  EXPECT_EQ(tree.min, min_kmer);
  EXPECT_EQ(tree.max, max_kmer);
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
  auto kmer2 = container::RI_Kmer{0};
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer2);
  EXPECT_NE(nullptr, tree.summary);
  EXPECT_EQ(tree.summary->size, 4);
  EXPECT_EQ(kmer1, tree.max);
  EXPECT_EQ(kmer2, tree.min);
  EXPECT_EQ(kmer2, tree.summary->max);
  EXPECT_EQ(kmer2, tree.summary->min);
  EXPECT_EQ(kmer1, tree.trees[0]->max);
  EXPECT_EQ(kmer1, tree.trees[0]->min);
}

TEST_F(VebTreeTest, FindOnOneInsert) {
  veb::Veb_tree tree{16};
  auto kmer1 = container::RI_Kmer{1};
  veb::insert(tree, kmer1);
  EXPECT_EQ(kmer1, veb::find(tree, kmer1));
}

TEST_F(VebTreeTest, FindOnOneInsertINT) {
  veb::Veb_tree tree{16};
  auto kmer1 = container::RI_Kmer{1};
  veb::insert(tree, kmer1);
  EXPECT_EQ(kmer1, veb::find(tree, kmer1.integer_rep));
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

TEST_F(VebTreeTest, FindOnTwoInsertINT) {
  veb::Veb_tree tree{16};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer0 = container::RI_Kmer{0};
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer0);
  EXPECT_EQ(kmer1, veb::find(tree, kmer1.integer_rep));
  EXPECT_EQ(kmer0, veb::find(tree, kmer0.integer_rep));
}

TEST_F(VebTreeTest, FindOnFiveInsert) {
  veb::Veb_tree tree{4};
  auto kmer0 = container::RI_Kmer{0};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer2 = container::RI_Kmer{2};
  auto kmer3 = container::RI_Kmer{3};
  veb::insert(tree, kmer0);
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer2);
  veb::insert(tree, kmer3);
  
  //std::cout << "Tree min " << tree.min.integer_rep << std::endl;
  //std::cout << "Tree max " << tree.max.integer_rep << std::endl;
  //std::cout << "Summary min " << tree.summary->min.integer_rep << std::endl;
  //std::cout << "Summary max " << tree.summary->max.integer_rep << std::endl;
  //std::cout << "Subtree0 min " << tree.trees[0]->max.integer_rep << std::endl;
  //std::cout << "Subtree0 max " << tree.trees[0]->max.integer_rep << std::endl;
  //std::cout << "Subtree1 min " << tree.trees[1]->max.integer_rep << std::endl;
  //std::cout << "Subtree1 max " << tree.trees[1]->max.integer_rep << std::endl;
  EXPECT_EQ(kmer0, veb::find(tree, kmer0));
  EXPECT_EQ(kmer1, veb::find(tree, kmer1));
  EXPECT_EQ(kmer2, veb::find(tree, kmer2));
  EXPECT_EQ(kmer3, veb::find(tree, kmer3));
}

TEST_F(VebTreeTest, FindOneSucc) {
  veb::Veb_tree tree{4};
  auto kmer0 = container::RI_Kmer{0};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer2 = container::RI_Kmer{2};
  auto kmer3 = container::RI_Kmer{3};
  veb::insert(tree, kmer0);
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer2);
  veb::insert(tree, kmer3);
  EXPECT_EQ(kmer1, veb::succ(tree, kmer0));
  EXPECT_EQ(kmer2, veb::succ(tree, kmer1));
  EXPECT_EQ(kmer3, veb::succ(tree, kmer2));
}

TEST_F(VebTreeTest, FindOneSuccINT) {
  veb::Veb_tree tree{4};
  auto kmer0 = container::RI_Kmer{0};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer2 = container::RI_Kmer{2};
  auto kmer3 = container::RI_Kmer{3};
  veb::insert(tree, kmer0);
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer2);
  veb::insert(tree, kmer3);
  EXPECT_EQ(kmer1, veb::succ(tree, kmer0.integer_rep));
  EXPECT_EQ(kmer2, veb::succ(tree, kmer1.integer_rep));
  EXPECT_EQ(kmer3, veb::succ(tree, kmer2.integer_rep));
}

TEST_F(VebTreeTest, Find4) {
  veb::Veb_tree tree{16};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer2 = container::RI_Kmer{2};
  auto kmer3 = container::RI_Kmer{3};
  auto kmer4 = container::RI_Kmer{4};
  auto kmer5 = container::RI_Kmer{5};
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer2);
  veb::insert(tree, kmer3);
  veb::insert(tree, kmer4);
  veb::insert(tree, kmer5);
  EXPECT_EQ(kmer2, veb::succ(tree, kmer1));
  EXPECT_EQ(kmer3, veb::succ(tree, kmer2));
  EXPECT_EQ(kmer4, veb::succ(tree, kmer3));
  EXPECT_EQ(kmer5, veb::succ(tree, kmer4));
}

TEST_F(VebTreeTest, FailFindINT) {
  veb::Veb_tree tree{4};
  auto kmer0 = container::RI_Kmer{0};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer2 = container::RI_Kmer{2};
  auto kmer3 = container::RI_Kmer{3};
  auto null_kmer = container::RI_Kmer(-1);
  veb::insert(tree, kmer0);
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer2);
  //veb::insert(tree, kmer3);
  EXPECT_EQ(kmer0, veb::find(tree, kmer0.integer_rep));
  EXPECT_EQ(kmer1, veb::find(tree, kmer1.integer_rep));
  EXPECT_EQ(null_kmer, veb::find(tree, kmer3.integer_rep));
}
