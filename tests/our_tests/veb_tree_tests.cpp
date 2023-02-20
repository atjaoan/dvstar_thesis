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
  veb::insert(tree, container::RI_Kmer{1});
  EXPECT_EQ(tree.min, container::RI_Kmer{1});
  EXPECT_EQ(tree.max, container::RI_Kmer{1});
  EXPECT_EQ(tree.size, 16);
  EXPECT_EQ(nullptr, tree.trees[0]);
  EXPECT_NE(nullptr, tree.summary);
}
TEST_F(VebTreeTest, SummaryTreeAfterTwoInsert) {
  veb::Veb_tree tree{16};
  veb::insert(tree, container::RI_Kmer{1});
  veb::insert(tree, container::RI_Kmer{2});
  EXPECT_NE(nullptr, tree.summary);
  EXPECT_EQ(tree.summary->size, 4);
  EXPECT_EQ(tree.summary->max, container::RI_Kmer{0});
  EXPECT_EQ(tree.summary->min, container::RI_Kmer{0});
}

TEST_F(VebTreeTest, SubTreeAfterTwoInsert) {
  veb::Veb_tree tree{16};
  veb::insert(tree, 1);
  veb::insert(tree, 0);
  EXPECT_NE(nullptr, tree.summary);
  EXPECT_EQ(container::RI_Kmer{0}, tree.summary->min);
  EXPECT_EQ(container::RI_Kmer{0}, tree.summary->max);
  EXPECT_EQ(tree.trees[0]->min, container::RI_Kmer{1});
  EXPECT_EQ(tree.trees[0]->max, container::RI_Kmer{1});
  EXPECT_EQ(tree.min, 0);
}

TEST_F(VebTreeTest, FindOnOneInsert) {
  veb::Veb_tree tree{16};
  veb::insert(tree, container::RI_Kmer{1});
  EXPECT_EQ(container::RI_Kmer{1}, veb::find(tree, container::RI_Kmer{1}));
}

TEST_F(VebTreeTest, FindOnTwoInsert) {
  veb::Veb_tree tree{16};
  veb::insert(tree, container::RI_Kmer{1});
  veb::insert(tree, container::RI_Kmer{0});
  EXPECT_EQ(container::RI_Kmer{1}, veb::find(tree, container::RI_Kmer{1}));
  EXPECT_EQ(container::RI_Kmer{0}, veb::find(tree, container::RI_Kmer{0}));
}

TEST_F(VebTreeTest, Insert16Find16) {
  std::vector<container::RI_Kmer> kmers {};
  for (size_t i = 0; i < 16; i++){
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