#pragma once 

#include <functional>
#include <filesystem>

#include "vlmc_from_kmers/kmer.hpp"

/*
  Stores VLMC (multiple k-mers) in a container. 
*/

using Kmer = vlmc::VLMCKmer; 

namespace container{
 
class VLMC_Container{

  public:
    VLMC_Container() = default;
    ~VLMC_Container() = default;
    VLMC_Container(const std::filesystem::path &path_to_bintree){}; 

    Kmer null_kmer{};
    virtual size_t size() const { return 0;};
    virtual void push(const Kmer &kmer){};
    virtual void for_each(const std::function<void(Kmer &kmer)> &){};
    virtual Kmer &get(const int i) { return null_kmer; };

};

class VLMC_vector : public VLMC_Container {

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

class VLMC_multi_vector : public VLMC_Container {

  private: 
    std::vector<Kmer> container{256}; 

  public: 
    VLMC_multi_vector() = default;
    ~VLMC_multi_vector() = default; 

    VLMC_multi_vector(const std::filesystem::path &path_to_bintree) {
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

    void push(const Kmer &kmer) override { 
      auto it_pos = container.begin() + get_index_rep(kmer);
      if(it_pos > container.end()){
        container.resize(get_index_rep(kmer) + 1);
      }
      container.insert(it_pos, kmer); 
      }

    void for_each(const std::function<void(Kmer &kmer)> &f) override {
      for (auto kmer : container){
        f(kmer);
      }
    }

    Kmer &get(const int i) override { return container[i]; }

    
  int get_index_rep(const Kmer &kmer) {
      int integer_value = 0;
      int offset = 1;
      for (int i = kmer.length - 1; i >= 0; i--) {
        auto kmer_2_bits = kmer.extract2bits(i) + 1;
        integer_value += (kmer_2_bits * offset);
        offset *= 4;
      }
      return integer_value;
    }
};

}

