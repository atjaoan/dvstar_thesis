#pragma once 

#include <functional>
#include <memory>
#include <limits.h>

#include <iostream>

#include "read_in_kmer.hpp"

/*
  Stores VLMC (multiple k-mers) in a van Emde Boas tree. 
*/


namespace veb{
 
struct Veb_tree{
  container::RI_Kmer min;
  container::RI_Kmer max;
  container::RI_Kmer null_kmer = container::RI_Kmer(-1);
  size_t size = 0;
  size_t nr_subtrees = 0;
  bool is_empty = true;
  std::shared_ptr<Veb_tree> summary = nullptr;
  std::vector<std::shared_ptr<Veb_tree>> trees{0, nullptr};

  Veb_tree() = default;
  // srqt(nr_trees) = |U|
  Veb_tree(size_t u) : is_empty{true}, min{INT_MAX}, max{-1} {
    if(u <= 2){
      summary = nullptr;
      size = 0;
      trees = std::vector<std::shared_ptr<Veb_tree>>{0, nullptr};
    } else{
      this->nr_subtrees = std::ceil(std::sqrt(u));
      this->size = u;
      summary = std::make_shared<Veb_tree>(nr_subtrees);
      trees = std::vector<std::shared_ptr<Veb_tree>>{nr_subtrees, nullptr};
    }
  };

  ~Veb_tree() = default;

  size_t tree_group(size_t in){
    if(nr_subtrees < 2) return 0;
    return (in / nr_subtrees); }
  size_t tree_group_index(size_t in){ return (in % nr_subtrees); }
  //container::RI_Kmer tree_group(container::RI_Kmer in){ return (in / size); }
  //container::RI_Kmer tree_group_index(container::RI_Kmer in){ return (in % size); }
};

void insert(Veb_tree &t, container::RI_Kmer in) {
  if(t.is_empty){
    t.min = in;
    t.max = in;
    t.is_empty = false;
    return;
  }
  if(in < t.min){
    container::RI_Kmer temp = t.min;
    t.min = in;
    in = temp;
  }
  if(in >= t.max){
    t.max = in;
  }
  if(t.size <= 2){
    return;
  }
  size_t c = t.tree_group(in.integer_rep);
  size_t i = t.tree_group_index(in.integer_rep);

  if(t.trees[c] == nullptr){
    t.trees[c] = std::make_shared<Veb_tree>(t.nr_subtrees);
  }

  if(t.trees[c]->is_empty){
    in.integer_rep = c;
    insert(*t.summary, in);
  }

  in.integer_rep = i;
  insert(*t.trees[c], in);
  return;
  }

container::RI_Kmer find(Veb_tree &t, container::RI_Kmer in) {
  if(t.min == in) return t.min;
  if(t.max == in) return t.max;

  if(t.size < in.integer_rep) return t.null_kmer;
  
  size_t c = t.tree_group(in.integer_rep);
  size_t i = t.tree_group_index(in.integer_rep);

  if(t.trees[c] == nullptr) return t.null_kmer;

  in.integer_rep = i;
  auto ret_kmer = find(*t.trees[c], in);
  if(ret_kmer.integer_rep != in.integer_rep) return t.null_kmer;

  ret_kmer.integer_rep += c*t.nr_subtrees;
  return ret_kmer;
}

container::RI_Kmer find(Veb_tree &t, int in) {
  if(t.min.integer_rep == in) return t.min;
  if(t.max.integer_rep == in) return t.max;

  if(t.size < in) return t.null_kmer;

  size_t c = t.tree_group(in);
  size_t i = t.tree_group_index(in);

  if(t.trees[c] == nullptr) return t.null_kmer;

  auto ret_kmer = find(*t.trees[c], i);
  if(ret_kmer.integer_rep != in) return t.null_kmer;

  ret_kmer.integer_rep += c*t.nr_subtrees;
  return ret_kmer;
}

container::RI_Kmer succ(Veb_tree&t, int in){
  if(t.is_empty) {
    return t.null_kmer;
  }
  if(in < t.min.integer_rep) {
    return t.min;
  }
  if(t.size <= 2){
    return t.max;
  }

  size_t c = t.tree_group(in);
  size_t i = t.tree_group_index(in);
  if(t.trees[c] == nullptr){
    return t.max;
  }
  if(i < t.trees[c]->max.integer_rep){
    in = i;
    auto ret_kmer = succ(*t.trees[c], in);
    ret_kmer.integer_rep += c*t.nr_subtrees;
    return ret_kmer; 
  } else{
    in = c;
    auto c_prim = succ(*t.summary, in);
    if(c_prim.is_null) {
      return t.max;
    } else {
      auto ret_kmer = t.trees[c_prim.integer_rep]->min;
      ret_kmer.integer_rep += c_prim.integer_rep*t.nr_subtrees;
      return ret_kmer;
    }
  }
}

container::RI_Kmer succ(Veb_tree&t, container::RI_Kmer in){
  if(t.is_empty) {
    return t.null_kmer;
  }
  if(in < t.min) {
    return t.min;
  }
  if(t.size <= 2){
    return t.max;
  }

  size_t c = t.tree_group(in.integer_rep);
  size_t i = t.tree_group_index(in.integer_rep);
  if(t.trees[c] == nullptr){
    return t.max;
  }
  if(i < t.trees[c]->max.integer_rep){
    in.integer_rep = i;
    auto ret_kmer = succ(*t.trees[c], in);
    ret_kmer.integer_rep += c*t.nr_subtrees;
    return ret_kmer; 
  } else{
    in.integer_rep = c;
    auto c_prim = succ(*t.summary, in);
    if(c_prim.is_null) {
      return t.max;
    } else {
      auto ret_kmer = t.trees[c_prim.integer_rep]->min;
      ret_kmer.integer_rep += c_prim.integer_rep*t.nr_subtrees;
      return ret_kmer;
    }
  }
}
}