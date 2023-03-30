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
  std::vector<int> container{};
  int size;
	int height;

  Veb_array() = default;
  Veb_array(int elems) : container(elems), size(elems), height{(int)(ceil(log((float)elems)/log(2.0)))} { };
  ~Veb_array() = default;

inline int power_of_two(int exponent) { return 1 << exponent; }

void hyper_compute(int n, int& d, int& D, int& subtree_size, int& subtree_leaf_count){	
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

void insert(int n){
  container[get_index(n, height)] = n;
}

int veb_search(const int *tree, int size, int elt){ 
  int d, D, subtree_size, subtree_leaf_count;
	hyper_compute(size, d, D, subtree_size, subtree_leaf_count);

	if (size > 1) {
		// Recurse on top half of tree
		int subtree_index = veb_search(tree, D, elt);
		if (subtree_index < 0) return subtree_index;

		int offset = subtree_index * subtree_size + D;

		// If not in top half, use subtree index to find place in bottom half
		int bottom_subtree_index = veb_search(tree + offset, subtree_size, elt);
		return subtree_leaf_count*subtree_index + bottom_subtree_index;
	} else {
		int root = tree[0];

		if (elt == root) {
			std::cout << "Found " << elt << " at index " 
				  << (int)(tree - container.data()) << std::endl;
			return -1;
		}

		return (elt < root) ? 0 : 1;		
	}
}
};
}