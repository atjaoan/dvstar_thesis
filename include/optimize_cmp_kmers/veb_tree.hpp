#pragma once 

#include <functional>
#include <memory>
#include <limits.h>

#include <iostream>

#include "read_in_kmer.hpp"

/*
  Stores VLMC (multiple k-mers) in a container. 
*/


namespace veb{
 
struct Veb_tree{
  container::RI_Kmer min;
  container::RI_Kmer max;
  container::RI_Kmer null_kmer = container::RI_Kmer(-1);
  size_t size = 0;
  bool is_empty;
  std::unique_ptr<Veb_tree> summary;
  std::vector<std::shared_ptr<Veb_tree>> trees;

  Veb_tree() = default;
  // srqt(nr_trees) = |U|
  Veb_tree(size_t nr_trees) : is_empty{true}, min{INT_MAX}, max{INT_MIN} {
    if(nr_trees <= 2){
      summary = nullptr;
      size = 0;
      trees = std::vector<std::shared_ptr<Veb_tree>>{0, nullptr};
    } else{
      size_t nr_subtrees = std::ceil(std::sqrt(nr_trees));
      this->size = nr_trees;
      summary = std::make_unique<Veb_tree>(nr_subtrees);
      trees = std::vector<std::shared_ptr<Veb_tree>>{nr_subtrees, nullptr};
    }
  };

  ~Veb_tree() = default;

  size_t tree_group(size_t in){ return (in / size); }
  size_t tree_group_index(size_t in){ return (in % size); }
  //container::RI_Kmer tree_group(container::RI_Kmer in){ return (in / size); }
  //container::RI_Kmer tree_group_index(container::RI_Kmer in){ return (in % size); }
};

void insert(Veb_tree &t, container::RI_Kmer& in) {
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
  if (!(in < t.max)){
    container::RI_Kmer temp = t.max;
    t.max = in;
    in = temp;
  }
  size_t c = t.tree_group(in.integer_rep);
  size_t i = t.tree_group_index(in.integer_rep);
  if(t.trees[c] == nullptr){
    in.integer_rep = c;
    insert(*t.summary, in);
  }

  if(t.trees[c] == nullptr)
    t.trees[c] = std::make_shared<Veb_tree>(t.size);
  in.integer_rep = i;
  insert(*t.trees[c], in);
  return;
  }

container::RI_Kmer& find(Veb_tree &t, container::RI_Kmer in) {
  if(t.min == in) return t.min;
  if(t.max == in) return t.max;
  size_t c = t.tree_group(in.integer_rep);
  find(*t.trees[c], in);
}

container::RI_Kmer succ(Veb_tree&t, container::RI_Kmer in){
  std::cout << "1" << std::endl;
  if(t.is_empty) return t.null_kmer;
  std::cout << "2" << std::endl;
  if(in < t.min) return t.min;
  std::cout << "3" << std::endl;

  size_t c = t.tree_group(in.integer_rep);
  size_t i = t.tree_group_index(in.integer_rep);
  auto temp = container::RI_Kmer(i);
  std::cout << "before if" << std::endl; 
  if(t.trees[c] != nullptr && in < t.trees[c]->max){
    std::cout << "first if" << std::endl;
    return succ(*t.trees[c], temp); //c*std::ceil(std::sqrt(t.size)) + 
  } else{
    std::cout << "else" << std::endl;
    int c_prim = succ(*t.summary, c).integer_rep;
    if(c_prim == -1) {
      std::cout << "c_prim == -1" << std::endl; 
      return t.max;
    } else {
      std::cout << "last" << std::endl; 
      return t.trees[c_prim]->min; //c_prim*std::ceil(std::sqrt(t.size)) + 
    }
  }
}
}