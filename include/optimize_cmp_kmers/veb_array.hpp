#pragma once 

#include <functional>
#include <memory>
#include <limits.h>
#include <iostream>

#include "read_in_kmer.hpp"

/*
  Stores VLMC (multiple k-mers) in a van Emde Boas layout. 
*/

namespace veb {
struct Veb_array {
  alignas(64) std::vector<container::RI_Kmer> a;
	container::RI_Kmer null_kmer = container::RI_Kmer(-1);
  static const unsigned MAX_H = 32;
  int height;
  int n;
  typedef unsigned char h_type;
	struct dumdum {
		h_type h0;
    h_type h1;
    h_type dummy[2];
		int m0;
    int m1;
	};
	dumdum s[MAX_H+1];

  Veb_array() = default;
  ~Veb_array() = default;

  void sequencer(int h, dumdum *s, unsigned d) {
	  if (h == 0) return;
	  int h0 = h/2;
	  int h1 = h-h0-1;
	  sequencer(h0, s, d);
	  s[d+h0].h0 = h0;
	  s[d+h0].m0 = (2<<h0)-1;
	  s[d+h0].h1 = h1;
	  s[d+h0].m1 = (2<<h1)-1;
	  sequencer(h1, s, d+h0+1);
  }

  container::RI_Kmer* construct(container::RI_Kmer* a0, int* rtl, int path, unsigned d) {
    if (d > height || rtl[d] >= n) return a0;

	  // visit left child
	  path <<= 1;
	  rtl[d+1] = rtl[d-s[d].h0] + s[d].m0 + (path&s[d].m0)*(s[d].m1);
	  a0 = construct(a0, rtl, path, d+1);

	  a[rtl[d]] = *a0++;

	  // visit right child
	  path += 1;
	  rtl[d+1] = rtl[d-s[d].h0] + s[d].m0 + (path&s[d].m0)*(s[d].m1);
	  a0 = construct(a0, rtl, path, d+1);

	  return a0;
  }

	Veb_array(std::vector<container::RI_Kmer>& from_container){
    n = from_container.size();
    // find smallest h such that sum_i=0^h 2^h >= n
	  int m = 1;
	  for (height = 0; m < n; height++, m += 1<<height);

    dumdum q = {(h_type)height, 0, {0, 0}, (2<<height)-1, 1};
	  std::fill_n(s, MAX_H+1, q);
	  sequencer(height, s, 0);

    a.reserve(n);
    int rtl[MAX_H + 1];
    rtl[0] = 0;
    construct(from_container.data(), rtl, 0, 0);
  }

  int search(int x) {
    int rtl[MAX_H+1];
	  int j = n;
	  int i = 0;
	  int p = 0;
	  for (int d = 0; i < n; d++) {
	  	rtl[d] = i;
	  	if (x < a[i].integer_rep) {
	  		p <<= 1;
	  		j = i;
	  	} else if (x > a[i].integer_rep) {
	  		p = (p << 1) + 1;
	  	} else {
	  		return i;
	  	}
	  	i = rtl[d-s[d].h0] + s[d].m0 + (p&s[d].m0)*(s[d].m1);
	  }
	  return j;
  }

  container::RI_Kmer& get_from_array(const int i_rep){
	  return a[search(i_rep)];
  }
  };
}

/*
namespace veb{

struct Veb_array{
  std::vector<container::RI_Kmer> container{};
	container::RI_Kmer null_kmer = container::RI_Kmer(-1);
  int size;
	int height;

  Veb_array() = default;
  Veb_array(int elems) : container(elems, null_kmer), size(elems), height{(int)(ceil(log((float)elems)/log(2.0)))} { };
  ~Veb_array() = default;

inline int power_of_two(int exponent) { return 1 << exponent; }

void get_group_index_size_count(int n, int& d, int& D, int& subtree_size, int& subtree_leaf_count){	
	int h = (int)ceil(log((float)n+1)/log(2.0f));
	d = h / 2;
	D = power_of_two(d) - 1;
	int delta = h - d;
	subtree_size = power_of_two(delta) - 1;
	subtree_leaf_count = power_of_two(h - 1);
}

int tree_size(int depth) { return power_of_two(depth) - 1; }

int get_index(int n, int height){ 
  if (height <= 1) return n;

	int bottom_half = height / 2;
	int top_half = height - bottom_half;

	int top_n = n >> bottom_half;
	int bottom_n = n & (power_of_two(bottom_half) - 1);

	if (bottom_n == 0) return get_index(top_n, top_half);

	int top_size = power_of_two(top_half) - 1;      
	int subtree_size = power_of_two(bottom_half) - 1;

	int top_address = top_n * tree_size(bottom_half) + tree_size(top_half);
	int bot_address = get_index(bottom_n, bottom_half);

	return top_address + bot_address;
}

container::RI_Kmer& retrive_on_index(const int i){
	return container[i];
}

void insert(const container::RI_Kmer kmer){
	int idx = get_index(kmer.integer_rep, height);
	if(idx > container.size()) {
		//std::cout << idx << " with size " << container.size() << "\n";
		return;
	}
  container[idx] = kmer;
}

container::RI_Kmer& get_elem(int n){
	int idx = get_index(n, height);
	if(idx > container.size() - 1) return null_kmer;
	return container[get_index(n, height)];
}

void print_array(){
	for(int i = 0; i < container.size(); i++){
		int rep = container[i].integer_rep;
		if(rep > 0){
			std::cout << rep << " at index " << i << "\n";
		}
	}
}
};
}
*/