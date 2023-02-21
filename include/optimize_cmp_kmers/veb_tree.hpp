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
  container::RI_Kmer null_kmer = container::RI_Kmer{};
  size_t size = 0;
  bool is_empty;
  std::unique_ptr<Veb_tree> summary;
  std::vector<std::shared_ptr<Veb_tree>> trees;

  Veb_tree() = default;
  // srqt(nr_trees) = |U|
  Veb_tree(size_t u) : is_empty{true}, min{INT_MAX}, max{INT_MIN} {
    if(u <= 2){
      summary = nullptr;
      size = 0;
      trees = std::vector<std::shared_ptr<Veb_tree>>{0, nullptr};
    } else{
      size_t nr_subtrees = std::ceil(std::sqrt(u));
      this->size = u;
      summary = std::make_unique<Veb_tree>(nr_subtrees);
      trees = std::vector<std::shared_ptr<Veb_tree>>{nr_subtrees, nullptr};
    }
  };

  ~Veb_tree() = default;

  size_t tree_group(size_t in){ return (in / std::sqrt(size)); }
  size_t tree_group_index(size_t in){ 
    std::cout << in << " % " << (int)std::sqrt(size) << std::endl;
    return (in % (int)std::sqrt(size)); }
  //container::RI_Kmer tree_group(container::RI_Kmer in){ return (in / size); }
  //container::RI_Kmer tree_group_index(container::RI_Kmer in){ return (in % size); }
};

void insert(Veb_tree &t, container::RI_Kmer in) {
  std::cout << "insert" << std::endl;
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
    container::RI_Kmer temp = t.max;
    t.max = in;
    //std::cout << "max is " << t.max.integer_rep << std::endl;
    in = temp;
  }
  if(t.size <= 2) return;
  //std::cout << in.integer_rep << " / " << t.size << std::endl;
  size_t c = t.tree_group(in.integer_rep);
  size_t i = t.tree_group_index(in.integer_rep);
  if(t.trees[c] == nullptr){
    t.trees[c] = std::make_shared<Veb_tree>((int)std::sqrt(t.size));
  }

  if(t.trees[c]->is_empty){
    in.integer_rep = c;
    insert(*t.summary, in);
  }
  in.integer_rep = i;
  insert(*t.trees[c], in);
  return;
  }

container::RI_Kmer& find(Veb_tree &t, container::RI_Kmer in) {
  if(t.min == in) return t.min;
  if(t.max == in) return t.max;
  size_t c = t.tree_group(in.integer_rep);
  in.integer_rep = c;
  find(*t.trees[c], in);
}

container::RI_Kmer& succ(Veb_tree&t, container::RI_Kmer in){
  //std::cout << "integer rep : " <<  in.integer_rep << std::endl;
  std::cout << "1" << std::endl;
  if(t.is_empty) {
    std::cout << "is empty " << std::endl;
    return t.null_kmer;
  }
  //std::cout << "2" << std::endl;
  if(in < t.min) return t.min;
  //std::cout << "3" << std::endl;
  if(t.size < 2)
    return t.max;

  size_t c = t.tree_group(in.integer_rep);
  size_t i = t.tree_group_index(in.integer_rep);
  if(t.trees[c] == nullptr){
    std::cout << "null -> tmax = " << t.max.integer_rep << std::endl;
    return t.max;
  }
  //std::cout << "tree max : " << t.trees[c]->max.integer_rep << std::endl;
  //std::cout << "i : " << i << std::endl;
  if(i < t.trees[c]->max.integer_rep){
    //std::cout << "first if" << std::endl;
    in.integer_rep = i;
    auto ret_kmer = succ(*t.trees[c], in);
    ret_kmer.integer_rep += c*(std::ceil)(std::sqrt(t.size));
    //std::cout << ret_kmer.integer_rep << std::endl;
    return ret_kmer; //c*std::ceil(std::sqrt(t.size)) + 
  } else{
    in.integer_rep = c;
    std::cout << "beefore " << std::endl;
    if(t.summary == nullptr)
      return t.max;
    auto c_prim = succ(*t.summary, in);
    std::cout << " after " << std::endl; 
    if(c_prim.is_null) {
      std::cout << " cprim " << std::endl; 
      return t.max;
    } else {
      std::cout << " hello " << std::endl;
      auto ret_kmer = t.trees[c_prim.integer_rep]->min;
      //std::cout << " 2 " << std::endl;
      ret_kmer.integer_rep += c_prim.integer_rep*std::ceil(std::sqrt(t.size));
      return ret_kmer; //c_prim*std::ceil(std::sqrt(t.size)) + 
    }
  }
}
}