#pragma once 

#include <iostream>

#include "read_in_kmer.hpp"
#include "global_aliases.hpp"

namespace b_tree{

class TreeNode {
  container::RI_Kmer *keys;
  int t;
  TreeNode **C;
  int n;
  bool leaf;
  
  public:

    TreeNode(int temp, bool bool_leaf);

    void insertNonFull(container::RI_Kmer k);
    void splitChild(int i, TreeNode *y);
    void traverse();

    container::RI_Kmer search(const int i_rep);

    void for_each(const std::function<void(const container::RI_Kmer &kmer)> &f);

    void second_pass(const eigenx_t &cached_context, 
              const size_t background_order, const size_t offset_to_remove);

    friend class BTree;
};

class BTree {
  TreeNode *root;
  int t;
  
  public:
    BTree(int temp) {
      root = NULL;
      t = temp;
    }

    void traverse() {
      if (root != NULL)
        root->traverse();
    }

    container::RI_Kmer search(const int i_rep) {
      if (root != NULL) {
        return root->search(i_rep);
      } 
      return Kmer{};
    }

    void for_each(const std::function<void(const container::RI_Kmer &kmer)> &f) {
      if (root != NULL) {
        root->for_each(f); 
      }
    }

    void second_pass(const eigenx_t &cached_context, 
              const size_t background_order, const size_t offset_to_remove) {
      if (root != NULL) {
        root->second_pass(cached_context, background_order, offset_to_remove);
      }
    }

    void insert(container::RI_Kmer k);
};

TreeNode::TreeNode(int t1, bool leaf1) {
  t = t1;
  leaf = leaf1;

  keys = new container::RI_Kmer[2 * t - 1];
  C = new TreeNode *[2 * t];

  n = 0;
}

void TreeNode::traverse() {
  int i;
  for (i = 0; i < n; i++) {
    if (leaf == false)
      C[i]->traverse();
    std::cout << " " << keys[i].integer_rep;
  }

  if (leaf == false)
    C[i]->traverse();
}

void TreeNode::for_each(const std::function<void(const container::RI_Kmer &kmer)> &f) {
  int i;
  for (i = 0; i < n; i++) {
    if (leaf == false)
      C[i]->for_each(f);
    f(keys[i]);
  }

  if (leaf == false)
    C[i]->for_each(f);
}

void TreeNode::second_pass(const eigenx_t &cached_context, 
              const size_t background_order, const size_t offset_to_remove) {
  int i;
  for (i = 0; i < n; i++) {
    if (leaf == false)
      C[i]->second_pass(cached_context, background_order, offset_to_remove);
    if(keys[i].is_null) {
      std::cout << "I am here" << std::endl;  
    } else {
      int background_idx = keys[i].background_order_index(keys[i].integer_rep, background_order);
      int offset = background_idx - offset_to_remove;
      keys[i].next_char_prob *= cached_context.row(offset).rsqrt();
    }
  }

  if (leaf == false)
    C[i]->second_pass(cached_context, background_order, offset_to_remove);
}

container::RI_Kmer TreeNode::search(const int i_rep) {
  int i = 0;
  while (i < n && keys[i].integer_rep < i_rep)
    i++;

  if (keys[i].integer_rep == i_rep)
    return keys[i];

  if (leaf == true)
    return container::RI_Kmer{};

  return C[i]->search(i_rep);
}

void BTree::insert(container::RI_Kmer k) {
  if (root == NULL) {
    root = new TreeNode(t, true);
    root->keys[0] = k;
    root->n = 1;
  } else {
    if (root->n == 2 * t - 1) {
      TreeNode *s = new TreeNode(t, false);

      s->C[0] = root;

      s->splitChild(0, root);

      int i = 0;
      if (s->keys[0] < k)
        i++;
      s->C[i]->insertNonFull(k);

      root = s;
    } else
      root->insertNonFull(k);
  }
}

void TreeNode::insertNonFull(container::RI_Kmer k) {
  int i = n - 1;

  if (leaf == true) {
    while (i >= 0 && k < keys[i]) {
      keys[i + 1] = keys[i];
      i--;
    }

    keys[i + 1] = k;
    n = n + 1;
  } else {
    while (i >= 0 && k < keys[i])
      i--;

    if (C[i + 1]->n == 2 * t - 1) {
      splitChild(i + 1, C[i + 1]);

      if (keys[i + 1] < k)
        i++;
    }
    C[i + 1]->insertNonFull(k);
  }
}

void TreeNode::splitChild(int i, TreeNode *y) {
  TreeNode *z = new TreeNode(y->t, y->leaf);
  z->n = t - 1;

  for (int j = 0; j < t - 1; j++)
    z->keys[j] = y->keys[j + t];

  if (y->leaf == false) {
    for (int j = 0; j < t; j++)
      z->C[j] = y->C[j + t];
  }

  y->n = t - 1;
  for (int j = n; j >= i + 1; j--)
    C[j + 1] = C[j];

  C[i + 1] = z;

  for (int j = n - 1; j >= i; j--)
    keys[j + 1] = keys[j];

  keys[i] = y->keys[t - 1];
  n = n + 1;
}
}