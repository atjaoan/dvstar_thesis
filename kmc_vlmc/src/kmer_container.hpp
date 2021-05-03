#pragma once

#include <execution>
#include <functional>
#include <memory>

#include <stxxl/sorter>

#include "kmer.hpp"

// template <int MAX_K>
// using kmer_sorter =
//     stxxl::sorter<VLMCKmer, ReverseKMerComparator<MAX_K>, 128 * 1024 * 1024>;

// template <int MAX_K> using kmer_sorter = stxxl::vector<VLMCKmer>;

template <class Comparator = ReverseKMerComparator<31>> class KmerContainer {
  static_assert(std::is_base_of<VirtualKMerComparator<31>, Comparator>::value,
                "Invalid comparator template");
  // using array_type = std::array<T, 3>;
  // using iterator = array_type::iterator;
  // using const_iterator = array_type::const_iterator;
public:
  KmerContainer() = default;
  ~KmerContainer() = default;

  virtual void push(VLMCKmer &kmer){};
  virtual void sort(){};
  virtual void for_each(const std::function<void(VLMCKmer &kmer)> &){};
};

template <class Comparator>
class InCoreKmerContainer : public KmerContainer<Comparator> {
  std::deque<VLMCKmer> container{};

public:
  InCoreKmerContainer() = default;
  ~InCoreKmerContainer() = default;

  void push(VLMCKmer &kmer) { container.push_back(kmer); }

  void sort() {
    std::sort(std::execution::par_unseq, container.begin(), container.end(),
              Comparator());
  };

  void for_each(const std::function<void(VLMCKmer &kmer)> &f) {
    std::for_each(container.begin(), container.end(), f);
  }
};

template <class Comparator>
class OutOfCoreKmerContainer : public KmerContainer<Comparator> {
  template <int MAX_K>
  using kmer_sorter =
      stxxl::sorter<VLMCKmer, VirtualKMerComparator<MAX_K>, 16 * 1024 * 1024>;
  kmer_sorter<31> sorter{Comparator(), 128 * 1024 * 1024};

public:
  OutOfCoreKmerContainer() = default;
  ~OutOfCoreKmerContainer() = default;

  void push(VLMCKmer &kmer) { sorter.push(kmer); }

  void sort() { sorter.sort(); };

  void for_each(std::function<void(VLMCKmer &kmer)> f) {
    sorter.rewind();

    VLMCKmer kmer{};

    while (!sorter.empty()) {
      kmer = *sorter;

      f(kmer);

      ++sorter;
    }
  }
};