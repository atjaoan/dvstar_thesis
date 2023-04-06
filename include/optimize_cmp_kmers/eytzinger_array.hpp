#pragma once 

#include <iostream>

#include "read_in_kmer.hpp"
#include <bits/stdc++.h>

namespace array {

struct Ey_array {
  //alignas(64) std::vector<container::RI_Kmer> a(size);
  //alignas(64) std::vector<container::RI_Kmer> b(size);
	container::RI_Kmer null_kmer = container::RI_Kmer(-1);
  int size = 1 << 20;
	int block_size = 2;
  alignas(64) container::RI_Kmer a[size], b[size+1];

  Ey_array() = default;
  Ey_array(std::vector<container::RI_Kmer>& from_container){
    std::sort(from_container.begin(), from_container.end());
    std::copy(from_container.begin(), from_container.end(), a);
  }
  ~Ey_array() = default;

int construct(int i = 0, int k = 1) {
    if (k <= n) {
        i = Ey_array::construct(i, 2 * k);
        b[k] = a[i++];
        i = Ey_array::construct(i, 2 * k + 1);
    }
    return i;
}

int search(int x) {
    int k = 1;
    while (k <= n) {
        __builtin_prefetch(b + k * block_size);
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