#pragma once 

#include <gtest/gtest.h>

#include <functional>
#include <filesystem>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "vlmc_from_kmers/kmer.hpp"
#include "read_in_kmer.hpp"
#include "global_aliases.hpp"
#include "b_tree.hpp"
#include "b_tree_layout.hpp"

using RI_Kmer = container::RI_Kmer; 

class BTreeTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path path_bintree_1{"../data/test_VLMCs/sequences_1.bintree"};
  std::filesystem::path path_bintree_2{"../data/test_VLMCs/sequences_2.bintree"};

  std::filesystem::path path_to_bintrees{"../data/test_VLMCs"};
};

b_tree::B_Tree createTree(std::filesystem::path path){
  std::ifstream ifs(path, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);
  std::vector<RI_Kmer> tmp{}; 
  Kmer kmer{};

  while (ifs.peek() != EOF){
    archive(kmer);
    RI_Kmer ri_kmer{kmer};
    tmp.push_back(ri_kmer);
  }
  ifs.close();
  std::sort(std::execution::seq, tmp.begin(), tmp.end());
  return b_tree::B_Tree(tmp);
}


TEST_F(BTreeTests, new_btree) {
}
