#pragma once 

#include <functional>
#include <filesystem>

#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_container.hpp"

/*
  Stores multiple VLMCs in a container. 
*/
namespace container {

template <typename VC> 
class Cluster_Container {

  private: 
    std::vector<VC> container{};

  public: 
    Cluster_Container() = default;
    ~Cluster_Container() = default; 

    Cluster_Container(const size_t i) : container(i) { }

    size_t size() const { return container.size(); }

    void push(const VC vlmc) { container.push_back(vlmc); }

    void for_each(const std::function<void(VC vlmc)> &f) {
      for (auto vlmc : container){
        f(vlmc);
      }
    }

    VLMC_Container &get(const int i) { return container[i]; }

    VC &operator[](size_t index){ return container[index]; }

    const VC &operator[](size_t index) const { return container[index]; }

};

struct Kmer_Pair {
  RI_Kmer kmer;
  size_t id; 

  Kmer_Pair(RI_Kmer kmer, size_t id){
    this->kmer = kmer; 
    this->id = id; 
  }

  Kmer_Pair() = default;

  ~Kmer_Pair() = default; 
};

class Kmer_Cluster {

  private: 
    std::unordered_map<int, Kmer_Pair> container{};

  public: 
    Kmer_Cluster() = default;
    ~Kmer_Cluster() = default; 

    void push(const Kmer_Pair kmer_pair) { container[kmer_pair.kmer.integer_rep] = kmer_pair; }

    std::vector<Kmer_Pair> &get(int bucket_num){
      std::vector<Kmer_Pair> returning_vector;
      for( auto it = container.begin(bucket_num); it!= container.end(bucket_num); ++it){
        returning_vector.push_back(it->second);
      }
      return returning_vector;
    }
};
}