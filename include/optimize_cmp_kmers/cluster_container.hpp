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

  using vlmc_container = container::VLMC_Container;

  private: 
    std::vector<VC> container{};

  public: 
    Cluster_Container() = default;
    ~Cluster_Container() = default; 

    size_t size() const { return container.size(); }

    void push(const VC vlmc) { container.push_back(vlmc); }

    void for_each(const std::function<void(VC vlmc)> &f) {
      for (auto vlmc : container){
        f(vlmc);
      }
    }

    vlmc_container &get(const int i) { return container[i]; }

};
}