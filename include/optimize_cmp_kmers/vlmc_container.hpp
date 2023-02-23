#pragma once 

#include <functional>
#include <filesystem>
#include <limits.h>
#include <exception>
#include "vlmc_from_kmers/kmer.hpp"
#include "optimize_cmp_kmers/read_in_kmer.hpp"
#include "b_tree.hpp"
#include "robin_hood.h"
#include "veb_tree.hpp"

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
    virtual RI_Kmer find(const int idx) { std::cout << "Bad place" << std::endl; return null_kmer; }
    virtual int get_max_kmer_index() const { return INT_MAX; }
    virtual int get_min_kmer_index() const { return 0; }

    virtual void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f){};
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
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f) override {
      for (size_t i = left_kmers.get_min_kmer_index() ; i <= left_kmers.get_max_kmer_index(); i++) {
        const RI_Kmer &left_kmer = left_kmers.get(i);
        if (left_kmer.is_null){
          continue; 
        }
        auto right_kmer = right_kmers.find(left_kmer.integer_rep);
        if (!right_kmer.is_null){
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
    int min_kmer_index = 0; // <- cant be set to MAX_INT when all kmers are inserted in order (0,1,2,3,4...)
    int container_size = -1; 
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
      if(index > container_size){
        container.resize((index + 1) * 2);
        max_kmer_index = index;
        container_size = (index + 1) * 2;
      } else if (index > max_kmer_index) {
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
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f) override {
      for (size_t i = left_kmers.get_min_kmer_index() ; i <= left_kmers.get_max_kmer_index(); i++) {
        const RI_Kmer &left_kmer = left_kmers.get(i);
        if (left_kmer.is_null){
          continue; 
        }
        auto right_kmer = right_kmers.find(left_kmer.integer_rep);
        if (!right_kmer.is_null){
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
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f) override {
      size_t left_i = 0;
      size_t right_i = 0;  
      size_t left_size = left_kmers.size();
      size_t right_size = right_kmers.size();
      while(left_i < left_size && right_i < right_size) {
        const RI_Kmer &left_kmer = left_kmers.get(left_i);
        const RI_Kmer &right_kmer = right_kmers.get(right_i);
        if (right_kmer == left_kmer) {
          f(left_kmer, right_kmer);
          left_i++;
          right_i++;
        } else if (left_kmer < right_kmer) {
          left_i++;
        } else {
          right_i++;
        }
      }
    }
};

/*
  Storing Kmers in a B-tree.
*/

class VLMC_B_tree : public VLMC_Container {

  private: 
    b_tree::BTree container{3};

  public: 
    VLMC_B_tree() = default;
    ~VLMC_B_tree() = default; 

    VLMC_B_tree(const std::filesystem::path &path_to_bintree) {
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

    size_t size() const override { return 0; }

    void push(const RI_Kmer &kmer) override { container.insert(kmer); }

    RI_Kmer &get(const int i) override { return null_kmer; }

    int get_max_kmer_index() const override { return -1; }
    int get_min_kmer_index() const override { return 0; }

    RI_Kmer find(const int i_rep) override {
      return container.search(i_rep); 
    }

    void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f) override {
      container.for_each([&](const RI_Kmer &left_v) {
        RI_Kmer right_v = right_kmers.find(left_v.integer_rep);
        if (!right_v.is_null){
          f(left_v, right_v);
        } 
      });
    }
};

/*
  Storing Kmers in a unordered map (HashMap).
*/
class VLMC_hashmap : public VLMC_Container {

  private: 
    robin_hood::unordered_map<int, RI_Kmer> container{};

  public: 
    VLMC_hashmap() = default;
    ~VLMC_hashmap() = default; 

    VLMC_hashmap(const std::filesystem::path &path_to_bintree) {
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

    void push(const RI_Kmer &kmer) override { container[kmer.integer_rep] = kmer; }

    RI_Kmer &get(const int i) override { return container[i]; }

    int get_max_kmer_index() const override { return container.size() - 1; }
    int get_min_kmer_index() const override { return 0; }

    RI_Kmer find(const int i_rep) override {
      auto res = container.find(i_rep);
      if (res != container.end()){
        return res->second; 
      }
      return null_kmer; 
    }

    void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &, const RI_Kmer &)> &f) override {
      for (auto &[i_rep, left_v] : container) {
        auto right_v = right_kmers.find(i_rep); 
        if (!right_v.is_null) {
          f(left_v, right_v); 
        } 
      }
    }
};


/*
  Storing Kmers in a vector where the Kmer string is used as index.
*/
class VLMC_Combo : public VLMC_Container {

  private: 
    const int max_idx = 4096; 
    std::array<RI_Kmer, 4096> container_ibv{};
    std::vector<RI_Kmer> container_sorted{};  
    int max_kmer_index = -1;
    int min_kmer_index = INT_MAX;
    RI_Kmer null_kmer {};

  public: 
    VLMC_Combo() = default;
    ~VLMC_Combo() = default; 

    VLMC_Combo(const std::filesystem::path &path_to_bintree) {
      std::ifstream ifs(path_to_bintree, std::ios::binary);
      cereal::BinaryInputArchive archive(ifs);

      Kmer input_kmer{};

      while (ifs.peek() != EOF){
        archive(input_kmer);
        RI_Kmer ri_kmer {input_kmer};
        push(ri_kmer);
      }
      ifs.close();
      std::sort(container_sorted.begin(), container_sorted.end());
    } 

    size_t size() const override { return container_sorted.size(); }

    void push(const RI_Kmer &kmer) override {  
      int index = kmer.integer_rep;
      if (index < max_idx){
        container_ibv[index] = kmer; 
      } else {
        container_sorted.push_back(kmer); 
      }
      }

    RI_Kmer &get(const int i) override { 
      return container_sorted[i]; 
    }

    int get_max_kmer_index() const override { return max_kmer_index; }
    int get_min_kmer_index() const override { return min_kmer_index; }

    RI_Kmer find(const int i_rep) override {
      if (i_rep < max_idx){
        return container_ibv[i_rep]; 
      } else {
        int L = 0;
        int R = size() - 1;
        // if (i_rep < R){
        //   R = i_rep;
        // }
        while (L <= R) {
            int m = (L + R) / 2;
            if (container_sorted[m].integer_rep < i_rep) {
              L = m + 1;
            } else if (container_sorted[m].integer_rep > i_rep) {
              R = m - 1;
            } else {
              return container_sorted[m];
            }
        }
      }
      return null_kmer;
    }

    void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f) override {
      for (size_t i = 0 ; i < max_idx; i++) {
        const RI_Kmer &left_kmer = container_ibv[i];
        if (left_kmer.is_null){
          continue; 
        }
        auto right_kmer = right_kmers.find(left_kmer.integer_rep);
        if (!right_kmer.is_null){
          f(left_kmer, right_kmer);
        }
      }
      size_t left_i = 0;
      size_t right_i = 0;  
      size_t left_size = left_kmers.size();
      size_t right_size = right_kmers.size();
      while(left_i < left_size && right_i < right_size) {
        const RI_Kmer &left_kmer = container_sorted[left_i];
        const RI_Kmer &right_kmer = right_kmers.get(right_i);
        if (right_kmer == left_kmer) {
          f(left_kmer, right_kmer);
          left_i++;
          right_i++;
        } else if (left_kmer < right_kmer) {
          left_i++;
        } else {
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

class VLMC_Veb : public VLMC_Container {

  private: 
    veb::Veb_tree veb = veb::Veb_tree(10000);

  public: 
    VLMC_Veb() = default;
    ~VLMC_Veb() = default; 

    VLMC_Veb(const std::filesystem::path &path_to_bintree) {
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

    size_t size() const override { return veb.size; }

    void push(const RI_Kmer &kmer) override { 
      if(veb.size < kmer.integer_rep){
        std::cout << "Too enourmous kmer : " << kmer.integer_rep << std::endl;
        return;
      }
      veb::insert(veb, kmer); 
    }

    RI_Kmer &get(const int i) override { return null_kmer; }

    int get_max_kmer_index() const override { return veb.size - 1; }
    int get_min_kmer_index() const override { return 0; }

    RI_Kmer find(const int i_rep) override { return veb::find(veb, i_rep); }

    void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f) override {
      RI_Kmer left_kmer = left_kmers.find(0);
      RI_Kmer right_kmer = right_kmers.find(0);

      //while(!left_kmer.is_null && !right_kmer.is_null)

      while(left_kmer.integer_rep != -1){
        // Check if right_kmers has this succeeding left_kmer
        right_kmer = right_kmers.find(left_kmer.integer_rep);
        if(right_kmer.integer_rep != -1){
          // if, apply f
          f(left_kmer, right_kmer);
        }
        // Iterating left
        left_kmer = veb::succ(veb, left_kmer);
      }
      // else continue
      /*
      for (size_t i = left_kmers.get_min_kmer_index() ; i <= left_kmers.get_max_kmer_index(); i++) {
        const RI_Kmer &left_kmer = left_kmers.get(i);
        if (left_kmer.is_null){
          continue; 
        }
        auto right_kmer = right_kmers.find(left_kmer.integer_rep);
        if (!right_kmer.is_null){
          f(left_kmer, right_kmer);
        } 
      }
      */
    }
};
}