#pragma once 

#include <functional>
#include <filesystem>

#include "vlmc_from_kmers/kmer.hpp"

/*
  Stores VLMC (multiple k-mers) in a container. 
*/

using Kmer = vlmc::VLMCKmer; 

namespace container{
 
class VLMC_template{

  public:
    VLMC_template() = default;
    ~VLMC_template() = default;
    VLMC_template(const std::filesystem::path &path_to_bintree){}; 

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

    VLMC_vector(const std::filesystem::path &path_to_bintree) {
      std::ifstream ifs(path_to_bintree, std::ios::binary);
      cereal::BinaryInputArchive archive(ifs);

      Kmer kmer{};

      while (ifs.peek() != EOF){
        archive(kmer);
        push(kmer);
      }
      ifs.close();
    } 

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

