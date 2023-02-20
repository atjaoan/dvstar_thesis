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
};

TEST_F(VebTreeTest, EmptyConstructor) {
  veb::Veb_tree tree{};
  EXPECT_EQ(tree.get_min(), INT_MAX);
  EXPECT_EQ(tree.get_max(), INT_MIN);
  EXPECT_EQ(tree.get_size(), 0);
  EXPECT_EQ(tree.get_summary(), nullptr);
  EXPECT_EQ(tree.get_trees(), nullptr);
  EXPECT_TRUE(tree.get_is_empty());
}

TEST_F(VebTreeTest, SizedConstructor) {
  veb::Veb_tree tree{16};
  EXPECT_EQ(tree.get_min(), INT_MAX);
  EXPECT_EQ(tree.get_max(), INT_MIN);
  EXPECT_EQ(tree.get_size(), 16);
  EXPECT_EQ(tree.get_summary(), nullptr);
  EXPECT_EQ(tree.get_trees()->get_size(), 0);
  EXPECT_TRUE(tree.get_is_empty());
}

TEST_F(VebTreeTest, InsertOne) {
  veb::Veb_tree tree{16};
  veb::insert(tree, 1);
  EXPECT_EQ(tree.get_min(), 1);
  EXPECT_EQ(tree.get_max(), 1);
  EXPECT_EQ(tree.get_size(), 16);
  EXPECT_EQ(tree.get_trees()->get_size(), 0);
}

TEST_F(VebTreeTest, SummaryTreeAfterTwoInsert) {
  veb::Veb_tree tree{16};
  veb::insert(tree, 1);
  veb::insert(tree, 2);
  EXPECT_NE(nullptr, tree.summary);
  EXPECT_EQ(tree.summary->size, 4);
  EXPECT_EQ(tree.summary->max, 0);
  EXPECT_EQ(tree.summary->min, 0);
}

TEST_F(VebTreeTest, SubTreeAfterTwoInsert) {
  veb::Veb_tree tree{16};
  veb::insert(tree, 1);
  veb::insert(tree, 0);
  EXPECT_NE(nullptr, tree.summary);
  EXPECT_EQ(0, tree.summary->min);
  EXPECT_EQ(0, tree.summary->max);
  EXPECT_EQ(tree.trees[0].min, 1);
  EXPECT_EQ(tree.trees[0].max, 1);
  EXPECT_EQ(tree.min, 0);
}