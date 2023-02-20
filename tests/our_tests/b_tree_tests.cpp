#pragma once 
#include <gtest/gtest.h>

#include <functional>
#include <filesystem>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "vlmc_from_kmers/kmer.hpp"
#include "read_in_kmer.hpp"
#include "b_tree.hpp"

class BTreeTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path path_bintree_1{"../data/test_VLMCs/sequences_1.bintree"};
  std::filesystem::path path_bintree_2{"../data/test_VLMCs/sequences_2.bintree"};
};

void fill_tree(b_tree::BTree &t, std::filesystem::path path){
  std::ifstream ifs(path, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);

  vlmc::VLMCKmer kmer{};

  while (ifs.peek() != EOF){
    archive(kmer);
    container::RI_Kmer ri_kmer{kmer};
    t.insert(ri_kmer);
  }
  ifs.close();
}

/*
void iterate_kmers(b_tree::TreeNode left_kmers, b_tree::TreeNode right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)>
        &f_not_shared) override {
  int left_i = 0; 
  int right_i = 0; 
  int left_n = left_kmers.n;
  int right_n = right_kmers.n; s 

  while (left_i < left_n && right_i < right_n){
    if (left_kmers == false)
  }
}

// void TreeNode::traverse() {
//   int i;
//   for (i = 0; i < n; i++) {
//     if (leaf == false)
//       C[i]->traverse();
//     std::cout << " " << keys[i].integer_rep;
//   }
// 
//   if (leaf == false)
//     C[i]->traverse();
// }
*/

TEST_F(BTreeTests, test) {
  b_tree::BTree t1(3);

  fill_tree(t1, path_bintree_1);

  b_tree::BTree t2(3);

  fill_tree(t2, path_bintree_2);

  std::cout << "The B-tree is: ";
  t1.traverse();
  std::cout << std::endl; 

  std::cout << "The B-tree is: ";
  t2.traverse();
  std::cout << std::endl; 
}
