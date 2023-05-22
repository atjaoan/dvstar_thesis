#pragma once 
#include <gtest/gtest.h>

#include <functional>
#include <filesystem>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <chrono>

#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_containers/veb_tree.hpp"
#include "vlmc_container.hpp"
#include "cluster_container.hpp"
#include "../read_helper.hpp"

void print_tree(veb::Veb_tree &t, int depth){
  std::cout << "Min at depth " << depth << " : " << t.min.integer_rep <<std::endl;
  std::cout << "Max at depth " << depth << " : " << t.max.integer_rep <<std::endl;
  if(t.summary != nullptr){
    std::cout << "Summary min at depth " << depth << " : " << t.summary->min.integer_rep <<std::endl;
    std::cout << "Summary max at depth " << depth << " : " << t.summary->max.integer_rep <<std::endl;
  } else {
    std::cout << "Null summary " << std::endl;
  }
  std::cout << "Size of t : " << t.size << std::endl;
  std::cout << "nr_subtrees : " << t.nr_subtrees << std::endl;
  if(t.trees[0] == nullptr){
    std::cout << "t.trees is null" << std::endl;
  //  return;
  }
  std::cout << "Subtrees " <<std::endl;
  for (int i = 0; i < t.trees.size(); i++)
  {
    if(t.trees[i] != nullptr){
      std::cout << "Subtree : " << i << std::endl;
      print_tree(*t.trees[i], depth +1);
    }
  }
}

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

TEST_F(VebTreeTest, FindOnFourInsert) {
  veb::Veb_tree tree{4};
  auto kmer0 = container::RI_Kmer{0};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer2 = container::RI_Kmer{2};
  auto kmer3 = container::RI_Kmer{3};
  veb::insert(tree, kmer0);
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer2);
  veb::insert(tree, kmer3);
  EXPECT_EQ(kmer0, veb::find(tree, kmer0));
  EXPECT_EQ(kmer1, veb::find(tree, kmer1));
  EXPECT_EQ(kmer2, veb::find(tree, kmer2));
  EXPECT_EQ(kmer3, veb::find(tree, kmer3));
}

TEST_F(VebTreeTest, FindThreeSucc) {
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
  veb::Veb_tree tree{10};
  auto kmer1 = container::RI_Kmer{0};
  auto kmer2 = container::RI_Kmer{1};
  auto kmer3 = container::RI_Kmer{2};
  auto kmer4 = container::RI_Kmer{3};
  auto kmer5 = container::RI_Kmer{4};
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
  EXPECT_EQ(kmer0, veb::find(tree, kmer0.integer_rep));
  EXPECT_EQ(kmer1, veb::find(tree, kmer1.integer_rep));
  EXPECT_EQ(null_kmer, veb::find(tree, kmer3.integer_rep));
}

TEST_F(VebTreeTest, Sized10) {
  veb::Veb_tree tree{10};
  auto kmer0 = container::RI_Kmer{0};
  auto kmer1 = container::RI_Kmer{1};
  auto kmer2 = container::RI_Kmer{2};
  auto kmer3 = container::RI_Kmer{3};
  auto null_kmer = container::RI_Kmer(-1);
  veb::insert(tree, kmer0);
  veb::insert(tree, kmer1);
  veb::insert(tree, kmer2);
  veb::insert(tree, kmer3);
  EXPECT_EQ(kmer1, veb::succ(tree, kmer0));
  EXPECT_EQ(kmer2, veb::succ(tree, kmer1));
  EXPECT_EQ(kmer3, veb::succ(tree, kmer2));
}

TEST_F(VebTreeTest, CheckMany) {
  for (size_t j = 10; j < 500; j++)
  {
  
  int items = j;
  veb::Veb_tree tree{items};

  for (size_t i = 0; i < items / 2; i++){
    auto kmer = container::RI_Kmer{i};
    veb::insert(tree, kmer);
  }

  int failed_finding_succ = 0;
  int found_nonexistent_succ = 0;

  for (size_t i = 0; i < items / 2 - 1; i++){
    container::RI_Kmer located = veb::succ(tree, i);
    if(located.integer_rep != i + 1){
      std::cout << "located : " << located.integer_rep << " should be : " << i+1 << std::endl;
      std::cout << "Failed at j = " << j << std::endl;
      std::cout << "With nr subtrees : " << std::ceil(std::sqrt(j)) << std::endl;
      failed_finding_succ++;
    }
  }
  for (size_t i = items / 2 + 1; i < items; i++){
    container::RI_Kmer located = veb::succ(tree, i);
    if(located.integer_rep != -1){
      found_nonexistent_succ++;
    }
  }
  EXPECT_EQ(0, failed_finding_succ);
  EXPECT_EQ(0, found_nonexistent_succ);
  }
}
