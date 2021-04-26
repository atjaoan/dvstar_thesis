#pragma once

#include <execution>
#include <functional>

#include "kmer.hpp"

// template <int MAX_K>
// using kmer_sorter =
//     stxxl::sorter<VLMCKmer, ReverseKMerComparator<MAX_K>, 128 * 1024 * 1024>;

// template <int MAX_K> using kmer_sorter = stxxl::vector<VLMCKmer>;

class KmerContainer {
  // using array_type = std::array<T, 3>;
  // using iterator = array_type::iterator;
  // using const_iterator = array_type::const_iterator;
public:
  virtual void push(VLMCKmer &kmer){};
  virtual void sort(){};
  virtual void for_each(std::function<void(VLMCKmer &kmer)>){};
};

class InCoreKmerContainer : public KmerContainer{
  std::deque<VLMCKmer> container{};

public:
  InCoreKmerContainer() = default;
  ~InCoreKmerContainer() = default;

  void push(VLMCKmer &kmer) { container.push_back(kmer); }

  void sort() {
    std::sort(std::execution::par_unseq, container.begin(), container.end(),
              ReverseKMerComparator<31>());
  };

  void for_each(std::function<void(VLMCKmer &kmer)> f) {
    std::for_each(container.begin(), container.end(), f);
  }
};

class OutOfCoreKmerContainer : public KmerContainer {
  template <int MAX_K>
  using kmer_sorter =
      stxxl::sorter<VLMCKmer, ReverseKMerComparator<MAX_K>, 16 * 1024 * 1024>;
  kmer_sorter<31> sorter{ReverseKMerComparator<31>(), 128 * 1024 * 1024};

public:
  OutOfCoreKmerContainer() = default;
  ~OutOfCoreKmerContainer() = default;

  void push(VLMCKmer &kmer) { sorter.push(kmer); }

  void sort() { sorter.sort(); };

  void for_each(std::function<void(VLMCKmer &kmer)> f) {
    VLMCKmer kmer{};

    while (!sorter.empty()) {
      kmer = *sorter;

      f(kmer);

      ++sorter;
    }
  }
};