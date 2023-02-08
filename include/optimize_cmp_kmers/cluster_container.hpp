#pragma once 

#include <functional>
#include <filesystem>

#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_template.hpp"

/*
  Stores multiple VLMCs in a container. 
*/

namespace cluster{
  
class Cluster_Container{

  using vlmc_template = container::VLMC_template;

  public:
    Cluster_Container() = default;
    ~Cluster_Container() = default;

    vlmc_template null_vlmc{};
    virtual size_t size() const { return 0; };
    virtual void push(const vlmc_template &vlmc) {};
    virtual void for_each(const std::function<void(vlmc_template &vlmc)> &){};
    virtual vlmc_template &get(const int i) { return null_vlmc; };

};

class Cluster_vector : public Cluster_Container {

  using vlmc_template = container::VLMC_template;

  private: 
    std::vector<vlmc_template> container{}; 

  public: 
    Cluster_vector() = default;
    ~Cluster_vector() = default; 

    size_t size() const override { return container.size(); }

    void push(const vlmc_template &vlmc) override { container.push_back(vlmc); }

    void for_each(const std::function<void(vlmc_template &vlmc)> &f) override {
      for (auto vlmc : container){
        f(vlmc);
      }
    }

    vlmc_template &get(const int i) override { return container[i]; }

};

}