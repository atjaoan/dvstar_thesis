#pragma once 

#include <iostream>

#include "read_in_kmer.hpp"
#include <bits/stdc++.h>

namespace array {

struct Ey_array {
	container::RI_Kmer null_kmer = container::RI_Kmer(-1);
  int size;
	const int block_size = 2; // 64 / sizeof(RI_Kmer)
  container::RI_Kmer *a;
  alignas(64) std::vector<container::RI_Kmer> b;

  Ey_array() = default;
  Ey_array(std::vector<container::RI_Kmer>& from_container){
    size = from_container.size();
    a = from_container.data();
    b.reserve(size + 1);
    b[0] = null_kmer;
    construct();
  }
  ~Ey_array() = default;

int construct(int i = 0, int k = 1) {
    if (k <= size) {
        i = Ey_array::construct(i, 2 * k);
        b[k] = a[i++];
        i = Ey_array::construct(i, 2 * k + 1);
    }
    return i;
}

int search(int x) {
    int k = 1;
    while (k <= size) {
        __builtin_prefetch(b.data() + k * block_size);
        k = 2 * k + (b[k] < x);
    }
    k >>= __builtin_ffs(~k);
    return k;
}

container::RI_Kmer& get_from_array(const int i_rep){
  return b[search(i_rep)];
}
};
}