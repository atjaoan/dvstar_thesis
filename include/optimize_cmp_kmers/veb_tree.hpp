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
  std::unique_ptr<Veb_tree[]> trees;

  Veb_tree() : trees{}, is_empty{true} {};
  // srqt(nr_trees) = |U|
  Veb_tree(size_t nr_trees) : 
    size{nr_trees}, trees{new Veb_tree[size]}, summary{}, is_empty{true} 
    {};
  ~Veb_tree() = default;
    
  size_t get_min() { return min; }
  size_t get_max() { return max; }
  size_t get_size() { return size; }

    //Temporary, to be removed
  Veb_tree* get_summary() { return summary.get(); }
  Veb_tree* get_trees() { return trees.get(); }
  bool get_is_empty() { return is_empty; }

  void remove(size_t elem) {}
  void pred(size_t elem) {}
  void succ(size_t elem) {}
    
  //Veb_tree(const std::filesystem::path &path_to_bintree){}; 

  //Kmer null_kmer{};
  //size_t size() const { return 0;};
  //void push(const Kmer &kmer){};
  //void for_each(const std::function<void(Kmer &kmer)> &){};
  //Kmer &get(const int i) { return null_kmer; };
  //std::tuple<std::reference_wrapper<Kmer>,bool> find(const Kmer &kmer) { return std::make_tuple(std::ref(null_kmer), false); }
  //int get_max_kmer_index() const { return INT_MAX; }
  //int get_min_kmer_index() const { return 0; }

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
  if(t.trees[c].is_empty){
    if(t.summary == nullptr) 
      t.summary = std::make_unique<Veb_tree>(std::sqrt(t.size));
    insert(*t.summary, c);
  }
  insert(t.trees[c], i);
  return;
  }

}
