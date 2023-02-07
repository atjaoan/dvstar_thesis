#pragma once 

#include <functional>
#include <filesystem>

#include "vlmc_from_kmers/kmer.hpp"

/*
  Stores multiple VLMCs in a container. 
*/

namespace cluster{
  
template<typename vlmc_c>
class Cluster_Container{

  public:
    Cluster_Container() = default;
    ~Cluster_Container() = default;
    // Cluster_Container(const std::filesystem::path &path_to_bintree); 

    vlmc_c null_vlmc{};
    virtual size_t size() const { return 0; };
    virtual void push(const vlmc_c &vlmc) {};
    virtual void for_each(const std::function<void(vlmc_c &vlmc)> &){};
    virtual vlmc_c &get(const int i) { return null_vlmc; };

};

template<typename vlmc_c> 
class Cluster_vector : public Cluster_Container<vlmc_c> {

  private: 
    std::vector<vlmc_c> container{}; 

  public: 
    Cluster_vector() = default;
    ~Cluster_vector() = default; 

    // Cluster_Container(const std::filesystem::path &path_to_bintree) 

    size_t size() const override { return container.size(); }

    void push(const vlmc_c &vlmc) override { container.push_back(vlmc); }

    void for_each(const std::function<void(vlmc_c &vlmc)> &f) override {
      for (auto vlmc : container){
        f(vlmc);
      }
    }

    vlmc_c &get(const int i) override { return container[i]; }

};

}