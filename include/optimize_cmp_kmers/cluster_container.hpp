#pragma once 

#include <functional>
#include <filesystem>

#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_container.hpp"

/*
  Stores multiple VLMCs in a container. 
*/

namespace container{
  
class Cluster_Container{

  using vlmc_container = container::VLMC_Container;

  public:
    Cluster_Container() = default;
    ~Cluster_Container() = default;

    vlmc_container null_vlmc{};
    virtual size_t size() const { return 0; };
    virtual void push(const vlmc_container &vlmc) {};
    virtual void for_each(const std::function<void(vlmc_container &vlmc)> &){};
    virtual vlmc_container &get(const int i) { return null_vlmc; };

};

class Cluster_vector : public Cluster_Container {

  using vlmc_container = container::VLMC_Container;

  private: 
    std::vector<vlmc_container> container{}; 

  public: 
    Cluster_vector() = default;
    ~Cluster_vector() = default; 

    size_t size() const override { return container.size(); }

    void push(const vlmc_container &vlmc) override { container.push_back(vlmc); }

    void for_each(const std::function<void(vlmc_container &vlmc)> &f) override {
      for (auto vlmc : container){
        f(vlmc);
      }
    }

    vlmc_container &get(const int i) override { return container[i]; }

};

}