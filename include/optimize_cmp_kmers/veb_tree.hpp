#pragma once 

#include <functional>
#include <memory>
#include <limits.h>

#include <iostream>

#include "vlmc_from_kmers/kmer.hpp"

/*
  Stores VLMC (multiple k-mers) in a container. 
*/

using Kmer = vlmc::VLMCKmer; 

namespace veb{
 
struct Veb_tree{
  size_t min = INT_MAX;;
  size_t max = INT_MIN;
  size_t size = 0;
  bool is_empty;
  std::unique_ptr<Veb_tree> summary;
  std::vector<std::shared_ptr<Veb_tree>> trees;

  Veb_tree() = default;
  // srqt(nr_trees) = |U|
  Veb_tree(size_t nr_trees) : is_empty{true} {
    if(nr_trees <= 2){
      summary = nullptr;
      size = 0;
      trees = std::vector<std::shared_ptr<Veb_tree>>{0, nullptr};
    } else{
      size_t nr_subtrees = std::ceil(std::sqrt(nr_trees));
      size = nr_subtrees;
      summary = std::make_unique<Veb_tree>(nr_subtrees);
      trees = std::vector<std::shared_ptr<Veb_tree>>{nr_subtrees, nullptr};
    }
  };

  ~Veb_tree() = default;

  size_t tree_group(size_t in){ return (in / size); }
  size_t tree_group_index(size_t in){ return (in % size); }
};

void insert(Veb_tree &t, size_t in) {
  if(t.is_empty){
    t.min = in;
    t.max = in;
    t.is_empty = false;
    return;
  }
  if(in < t.min){
    int temp = t.min;
    t.min = in;
    in = temp;
  } 
  if (in >= t.max){
    int temp = t.max;
    t.max = in;
    in = temp;
  }
  size_t c = t.tree_group(in);
  size_t i = t.tree_group_index(in);
  if(t.trees[c] == nullptr){
    insert(*t.summary, c);
  }

  if(t.trees[c] == nullptr)
    t.trees[c] = std::make_shared<Veb_tree>(t.size);
  insert(*t.trees[c], i);
  return;
  }

size_t find(Veb_tree &t, size_t in) {
  if(t.min == in) return t.min;
  if(t.max == in) return t.max;
  size_t c = t.tree_group(in);
  find(*t.trees[c], in);
}
}