#pragma once 

#include <functional>
#include <filesystem>

#include "vlmc_from_kmers/kmer.hpp"

/*
  Stores multiple VLMCs in a container. 
*/

using Kmer = vlmc::VLMCKmer; 

namespace container{
 
class VLMC_template{

  public:
    VLMC_template() = default;
    ~VLMC_template() = default;
    VLMC_template(const std::filesystem::path &directory); 

    Kmer null_kmer{};
    virtual size_t size() const { return 0; };
    virtual void push(const Kmer &kmer){};
    virtual void for_each(const std::function<void(Kmer &kmer)> &){};
    virtual Kmer &get(const int i) { return null_kmer; };

};

class VLMC_vector : public VLMC_template {

  private: 
    std::vector<Kmer> container{}; 

  public: 
    VLMC_vector() = default;
    ~VLMC_vector() = default; 
    VLMC_vector(const std::filesystem::path &directory); 

    size_t size() const override { return container.size(); }

    void push(const Kmer &kmer) override { container.push_back(kmer); }

    void for_each(const std::function<void(Kmer &kmer)> &f) override {
      for (auto kmer : container){
        f(kmer);
      }
    }

    Kmer &get(const int i) override { return container[i]; }

};

}

