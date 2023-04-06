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