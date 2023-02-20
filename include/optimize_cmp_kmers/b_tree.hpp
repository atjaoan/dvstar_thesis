// Searching a key on a B-tree in C++

#include <iostream>

#include "read_in_kmer.hpp"

namespace b_tree{

using Kmer = container::RI_Kmer;

class TreeNode {
  Kmer *keys;
  int t;
  TreeNode **C;
  int n;
  bool leaf;
  
  public:

    TreeNode(int temp, bool bool_leaf);

    void insertNonFull(Kmer k);
    void splitChild(int i, TreeNode *y);
    void traverse();

    Kmer search(const int i_rep);

    void for_each(const std::function<void(const Kmer &kmer)> &f);

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

    Kmer search(const int i_rep) {
      if (root != NULL) {
        return root->search(i_rep);
      } 
      return Kmer{};
    }

    void for_each(const std::function<void(const Kmer &kmer)> &f) {
      if (root != NULL) {
        root->for_each(f); 
      }
    }

    void insert(Kmer k);
};

TreeNode::TreeNode(int t1, bool leaf1) {
  t = t1;
  leaf = leaf1;

  keys = new Kmer[2 * t - 1];
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

void TreeNode::for_each(const std::function<void(const Kmer &kmer)> &f) {
  int i;
  for (i = 0; i < n; i++) {
    if (leaf == false)
      C[i]->for_each(f);
    f(keys[i]);
  }

  if (leaf == false)
    C[i]->for_each(f);
}

Kmer TreeNode::search(const int i_rep) {
  int i = 0;
  while (i < n && keys[i].integer_rep < i_rep)
    i++;

  if (keys[i].integer_rep == i_rep)
    return keys[i];

  if (leaf == true)
    return Kmer{};

  return C[i]->search(i_rep);
}

void BTree::insert(Kmer k) {
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

void TreeNode::insertNonFull(Kmer k) {
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