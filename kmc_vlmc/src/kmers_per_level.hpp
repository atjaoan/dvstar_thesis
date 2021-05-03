#pragma once

#include <vector>

template <class T> struct KMersPerLevel {
  int size_;
  int n_levels;
  std::vector<std::vector<T>> levels;

  KMersPerLevel(int size__, int n_levels_) : size_(size__), n_levels(n_levels) {
    levels = std::vector<std::vector<T>>(n_levels_, std::vector<T>(size__));
  };

  std::vector<T> &operator[](size_t pos) { return levels[pos]; }

  std::vector<T> operator[](size_t pos) const { return levels[pos]; }

  void reset_depth(int depth) {
    for (int i = 0; i < 4; i++) {
      levels[depth][i] = T{};
    }
  }

  size_t size() const { return size_; }
};