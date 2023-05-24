#pragma once 

#include <iostream>

#include "read_in_kmer.hpp"
#include <bits/stdc++.h>

namespace array {

struct Ey_array {
	kmers::RI_Kmer null_kmer = kmers::RI_Kmer(-1);
  int size;
	static const int block_size = 2; // = 64 / sizeof(RI_Kmer)
  kmers::RI_Kmer *kmer_from;
  alignas(64) std::vector<kmers::RI_Kmer> ey_sorted_kmers;

  Ey_array() = default;
  Ey_array(std::vector<kmers::RI_Kmer>& from_container){
    size = from_container.size();
    kmer_from = from_container.data();
    ey_sorted_kmers.reserve(size + 1);
    ey_sorted_kmers[0] = null_kmer;
    construct();
  }
  ~Ey_array() = default;

int construct(int i = 0, int k = 1) {
    if (k <= size) {
        i = Ey_array::construct(i, 2 * k);
        ey_sorted_kmers[k] = kmer_from[i++];
        i = Ey_array::construct(i, 2 * k + 1);
    }
    return i;
}

int search(int x) {
    int k = 1;
    while (k <= size) {
        __builtin_prefetch(ey_sorted_kmers.data() + k * block_size);
        k = 2 * k + (ey_sorted_kmers[k] < x);
    }
    k >>= __builtin_ffs(~k);
    return k;
}

kmers::RI_Kmer& get_from_array(const int i_rep){
  return ey_sorted_kmers[search(i_rep)];
}
};
}

namespace integer_array {

struct Min_max_node {
  int block_start;
  int max; // min, max; 

  Min_max_node(int idx, int max){
    this->block_start = idx; 
    // this->min = min;
    this->max = max;
  }

  Min_max_node() = default;
  ~Min_max_node() = default; 
};

struct Ey_integer_array {
  Min_max_node null_node;
  int size;
  int skip_size;
	static const int block_size = 8; // = 64 / sizeof(RI_Kmer)
  kmers::RI_Kmer *kmer_from;
  alignas(64) std::vector<Min_max_node> ey_sorted_int;

  Ey_integer_array() = default;
  Ey_integer_array(std::vector<kmers::RI_Kmer>& from_container){
    skip_size = std::ceil(std::log(from_container.size()));
    size = from_container.size() / skip_size;
    null_node = Min_max_node(from_container.size(), -1); 
    kmer_from = from_container.data();
    ey_sorted_int.reserve(size + 1);
    ey_sorted_int[0] = null_node;
    construct();
  }
  ~Ey_integer_array() = default;

int construct(int i = 0, int k = 1) {
    if (k <= size) {
        i = construct(i, 2 * k);
        ey_sorted_int[k] = Min_max_node(i, kmer_from[i+skip_size - 1].integer_rep);
        i += skip_size; 
        i = construct(i, 2 * k + 1);
    }
    return i;
}

int search(int x) {
    int k = 1;
    while (k <= size) {
        __builtin_prefetch(ey_sorted_int.data() + k * block_size);
        k = 2 * k + (ey_sorted_int[k].max < x);
    }
    k >>= __builtin_ffs(~k);
    return k;
}

Min_max_node get_from_array(const int i_rep){
  auto idx = search(i_rep); 
  if (ey_sorted_int[idx].max < i_rep) {
    while (true) {
      idx++; 
      if (idx >= size) {
        return null_node; 
      } else if (ey_sorted_int[idx].max > i_rep) {
        break; 
      }
    }
  } 
  return ey_sorted_int[idx]; 
}
};  
}