#pragma once 

#include <functional>
#include <filesystem>
#include <limits.h>

#include "vlmc_from_kmers/kmer.hpp"
#include "optimize_cmp_kmers/read_in_kmer.hpp"

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

    RI_Kmer null_kmer{};
    virtual size_t size() const { return 0;};
    virtual void push(const RI_Kmer &kmer){};
    virtual RI_Kmer &get(const int i) { std::cout << "Hello from bad place" << std::endl; return null_kmer; };
    virtual RI_Kmer find(const int idx) { 
      std::cout << "Bad place" << std::endl; 
      return null_kmer; 
    }
    virtual int get_max_kmer_index() const { return INT_MAX; }
    virtual int get_min_kmer_index() const { return 0; }

    virtual void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)>
        &f_not_shared){};
};

/*
  Storing Kmers in a unsorted vector.
*/
class VLMC_vector : public VLMC_Container {

  private: 
    std::vector<RI_Kmer> container{}; 

  public: 
    VLMC_vector() = default;
    ~VLMC_vector() = default; 

    VLMC_vector(const std::filesystem::path &path_to_bintree) {
      std::ifstream ifs(path_to_bintree, std::ios::binary);
      cereal::BinaryInputArchive archive(ifs);

      Kmer kmer{};

      while (ifs.peek() != EOF){
        archive(kmer);
        RI_Kmer ri_kmer{kmer};
        push(ri_kmer);
      }
      ifs.close();
    } 

    size_t size() const override { return container.size(); }

    void push(const RI_Kmer &kmer) override { container.push_back(kmer); }

    RI_Kmer &get(const int i) override { return container[i]; }

    int get_max_kmer_index() const override { return container.size() - 1; }
    int get_min_kmer_index() const override { return 0; }

    RI_Kmer find(const int i_rep) override {
      for (size_t i = 0; i < container.size(); i++){
        if (container[i].integer_rep==i_rep) {
          return container[i]; 
        }
      }
      return null_kmer; 
    }

    void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)>
        &f_not_shared) override {
      for (size_t i = left_kmers.get_min_kmer_index() ; i <= left_kmers.get_max_kmer_index(); i++) {
        const RI_Kmer &left_kmer = left_kmers.get(i);
        if (left_kmer.is_null){
          continue; 
        }
        auto right_kmer = right_kmers.find(left_kmer.integer_rep);
        if (right_kmer.is_null){
          f_not_shared(left_kmer, right_kmer);
        } else {
          f(left_kmer, right_kmer);
        }
      }
    }
};

/*
  Storing Kmers in a vector where the Kmer string is used as index.
*/
class Index_by_value : public VLMC_Container {

  private: 
    std::vector<RI_Kmer> container{}; 
    int c_size = 0;
    int max_kmer_index = -1;
    int min_kmer_index = INT_MAX;
    RI_Kmer null_kmer {};

  public: 
    Index_by_value() = default;
    ~Index_by_value() = default; 

    Index_by_value(const std::filesystem::path &path_to_bintree) {
      std::ifstream ifs(path_to_bintree, std::ios::binary);
      cereal::BinaryInputArchive archive(ifs);

      Kmer input_kmer{};

      while (ifs.peek() != EOF){
        archive(input_kmer);
        RI_Kmer ri_kmer {input_kmer};
        push(ri_kmer);
      }
      ifs.close();
      
    } 

    size_t size() const override { return c_size; }

    void push(const RI_Kmer &kmer) override {  
      int index = kmer.integer_rep;
      if(index > max_kmer_index){
        container.resize(index + 10);
        max_kmer_index = index;
      } else if (index < min_kmer_index){
        min_kmer_index = index;
      }
      //Must be done after resize (resize invalidades all iterators)
      // std::cout << "Index " << index << " maxKmerIndex = " << max_kmer_index << std::endl;
      container[index] = kmer;
      c_size++;
      }

    RI_Kmer &get(const int i) override { return container[i]; }

    int get_max_kmer_index() const override { return max_kmer_index; }
    int get_min_kmer_index() const override { return min_kmer_index; }

    RI_Kmer find(const int i_rep) override {
      if (i_rep <= max_kmer_index){
        return container[i_rep]; 
      }
      return null_kmer;
    }

    void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)>
        &f_not_shared) override {
      for (size_t i = left_kmers.get_min_kmer_index() ; i <= left_kmers.get_max_kmer_index(); i++) {
        const RI_Kmer &left_kmer = left_kmers.get(i);
        if (left_kmer.is_null){
          continue; 
        }
        auto right_kmer = right_kmers.find(left_kmer.integer_rep);
        if (right_kmer.is_null){
          f_not_shared(left_kmer, right_kmer);
        } else {
          f(left_kmer, right_kmer);
        }
      }
    }
};

/*
  Storing Kmers in a sorted vector.
*/
class VLMC_sorted_vector : public VLMC_Container {

  private: 
    std::vector<RI_Kmer> container{}; 

  public: 
    VLMC_sorted_vector() = default;
    ~VLMC_sorted_vector() = default; 

    VLMC_sorted_vector(const std::filesystem::path &path_to_bintree) {
      std::ifstream ifs(path_to_bintree, std::ios::binary);
      cereal::BinaryInputArchive archive(ifs);

      Kmer kmer{};

      while (ifs.peek() != EOF){
        archive(kmer);
        RI_Kmer ri_kmer{kmer};
        push(ri_kmer);
      }
      ifs.close();

      std::sort(container.begin(), container.end());
    } 

    size_t size() const override { return container.size(); }

    void push(const RI_Kmer &kmer) override { container.push_back(kmer); }

    RI_Kmer &get(const int i) override { return container[i]; }

    int get_max_kmer_index() const override { return container.size() - 1; }
    int get_min_kmer_index() const override { return 0; }

    RI_Kmer find(const int i_rep) override {
      int L = 0;
      int R = size() - 1;
      if (i_rep < R){
        R = i_rep;
      }
      while (L <= R) {
          int m = (L + R) / 2;
          if (container[m].integer_rep < i_rep) {
            L = m + 1;
          } else if (container[m].integer_rep > i_rep) {
            R = m - 1;
          } else {
            return container[m];
          }
      }
      return null_kmer;
    }

    void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)>
        &f_not_shared) override {
      size_t left_i = 0;
      size_t right_i = 0;  
      while(left_i < left_kmers.size() && right_i < right_kmers.size()) {
        const RI_Kmer &left_kmer = left_kmers.get(left_i);
        const RI_Kmer &right_kmer = right_kmers.get(right_i);
        if (right_kmer == left_kmer) {
          f(left_kmer, right_kmer);
          left_i++;
          right_i++;
        } else if (left_kmer < right_kmer) {
          f_not_shared(left_kmer, right_kmer);
          left_i++;
        } else {
          f_not_shared(right_kmer, left_kmer);
          right_i++;
        }
      }
    }
};

/*
class VLMC_multi_vector : public VLMC_Container {

  private: 
    std::vector<Kmer> container{}; 
    int c_size = 0;
    int max_kmer_index = 0;
    int min_kmer_index = INT_MAX;

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

    size_t size() const override { return c_size; }

    void push(const Kmer &kmer) override { 
      int index = get_index_rep(kmer);
      if(index > max_kmer_index){
        container.resize(index + 10);
        max_kmer_index = index;
      } else if (index < min_kmer_index){
        min_kmer_index = index;
      }
      //Must be done after resize (resize invalidades all iterators)
      container[index] = kmer; 
      c_size = c_size + 1;
      }

    void for_each(const std::function<void(Kmer &kmer)> &f) override {
      for (auto kmer : container){
        f(kmer);
      }
    }

    Kmer &get(const int i) override { return container[i]; }

    int get_max_kmer_index() const override { return max_kmer_index; }
    int get_min_kmer_index() const override { return min_kmer_index; }

    std::tuple<std::reference_wrapper<Kmer>,bool> find(const Kmer &kmer) override {
      auto index = get_index_rep(kmer);
      if (index <= max_kmer_index){
        if (container[index]==kmer){
          return std::make_tuple(std::ref(container[index]), true);
        }
      }
      return std::make_tuple(std::ref(null_kmer), false); 
    }

    int get_index_rep(const RI_Kmer &kmer) {
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
*/

}

