#pragma once 

#include <functional>
#include <filesystem>
#include <limits.h>
#include <exception>
#include <algorithm>
#include <set>
#include <execution>
#include "vlmc_from_kmers/kmer.hpp"
#include "optimize_cmp_kmers/read_in_kmer.hpp"
#include "b_tree.hpp"
#include "robin_hood.h"
#include "unordered_dense.h"
#include "veb_tree.hpp"
#include "Eigen/Core"

/*
  Stores VLMC (multiple k-mers) in a container. 
*/
constexpr int misses_before_skip = 6;

using Kmer = vlmc::VLMCKmer; 

namespace container{
 
class VLMC_Container{

  public:
    VLMC_Container() = default;
    ~VLMC_Container() = default;
    VLMC_Container(const std::filesystem::path &path_to_bintree, const size_t background_order){}; 

    RI_Kmer null_kmer{};
    virtual size_t size() const { return 0;};
    virtual void push(const RI_Kmer &kmer){};
    virtual RI_Kmer &get(const int i) { std::cout << "Hello from bad place" << std::endl; return null_kmer; };
    virtual RI_Kmer find(const int idx) { std::cout << "Bad place" << std::endl; return null_kmer; }
    virtual int get_max_kmer_index() const { return INT_MAX; }
    virtual int get_min_kmer_index() const { return 0; }
    virtual std::vector<RI_Kmer>::iterator begin() { return std::vector<RI_Kmer>{}.begin(); };
    virtual std::vector<RI_Kmer>::iterator end() { return std::vector<RI_Kmer>{}.end(); };

    virtual void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f){};
};

int load_VLMCs_from_file(const std::filesystem::path &path_to_bintree, Eigen::ArrayX4f &cached_context, 
        const std::function<void(const RI_Kmer &kmer)> f, const size_t background_order = 0) {
  std::ifstream ifs(path_to_bintree, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);
  Kmer input_kmer{};

  int offset_to_remove = 0;
  for (int i = 0; i < background_order; i++){
    offset_to_remove += std::pow(4, i); 
  }

  while (ifs.peek() != EOF){
    archive(input_kmer);
    RI_Kmer ri_kmer{input_kmer};
    if(input_kmer.length <= background_order){
      if (input_kmer.length + 1 > background_order){
        int offset = ri_kmer.integer_rep - offset_to_remove; 
        cached_context.row(offset) = ri_kmer.next_char_prob;
      }
    } else {
      f(ri_kmer); 
    }
  }
  ifs.close();

  return offset_to_remove; 
}

/*
  Storing Kmers in a unsorted vector.
*/
class VLMC_vector : public VLMC_Container {

  private: 
    std::vector<RI_Kmer> container{}; 

  public: 
    VLMC_vector() = default;
    ~VLMC_vector() = default; 

    VLMC_vector(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      Eigen::ArrayX4f cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      // Second pass - Comment this and select in dvstar.hpp row 120 to 128 to use old or new implementation. 
      for (size_t i = 0; i <= get_max_kmer_index(); i++){
        RI_Kmer kmer = get(i);
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        get(i).next_char_prob *= cached_context.row(offset).rsqrt();
      }
    } 

    // RI_Kmer* begin(){
    //   return &container[0];
    // }
// 
    // RI_Kmer* end(){
    //   return begin() + container.size();
    // }

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

class VLMC_Indexing : public VLMC_Container { 

  private: 
    std::vector<RI_Kmer> container{}; 
    int c_size = 0;
    int max_kmer_index = 0;
    int min_kmer_index = 0; // <- cant be set to MAX_INT when all kmers are inserted in order (0,1,2,3,4...)
    int container_size = -1; 

  public: 
    VLMC_Indexing(const int initial_size = 50) : container(initial_size, null_kmer), container_size{initial_size} {}
    ~VLMC_Indexing() = default; 

    VLMC_Indexing(const std::filesystem::path &path_to_bintree, const size_t background_order = 0, const int initial_size = 50) 
      : container(initial_size, null_kmer), container_size{initial_size} {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      Eigen::ArrayX4f cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      // Second Pass
      for (size_t i = 0; i <= max_kmer_index; i++){
        RI_Kmer kmer = container[i];
        if(kmer.is_null) continue;
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        get(i).next_char_prob *= cached_context.row(offset).rsqrt();
      }
    } 

    size_t size() const override { return c_size; }

    void push(const RI_Kmer &kmer) override {  
      int index = kmer.integer_rep;
      if(index >= container_size){
        container.resize((index + 1) * 2, null_kmer);
        max_kmer_index = index;
        container_size = (index + 1) * 2;
      } else if (index > max_kmer_index) {
        max_kmer_index = index; 
      } else if (index < min_kmer_index){
        min_kmer_index = index;
      }
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

    VLMC_sorted_vector(const std::filesystem::path &path_to_bintree, const size_t background_order = 0, bool use_new = false) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      Eigen::ArrayX4f cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);
      
      std::sort(std::execution::seq, container.begin(), container.end());
      for (size_t i = 0; i < size(); i++){
        RI_Kmer kmer = get(i);
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        get(i).next_char_prob *= cached_context.row(offset).rsqrt();
      }
    } 

    size_t size() const override { return container.size(); }

    void push(const RI_Kmer &kmer) override { container.push_back(kmer); }

    std::vector<RI_Kmer>::iterator begin() override { return container.begin(); };
    std::vector<RI_Kmer>::iterator end() override { return container.end(); };

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
      auto right_it = right_kmers.begin();
      auto right_end = right_kmers.end();
      
      auto left_it = left_kmers.begin();
      auto left_end = left_kmers.end();
      while(left_it != left_end && right_it != right_end){
        if(*left_it == *right_it){
          f(*left_it, *right_it);
          ++left_it;
          ++right_it;
        } else if(*left_it < *right_it) {
          ++left_it;
        }
        else ++right_it;
      }
      /*
      for(auto &left_kmer : left_kmers){
        while((*right_it) < left_kmer){
          ++right_it;
          if(right_it == right_end) break;
        }
        if(left_kmer == (*right_it)){
          f(left_kmer, *right_it);
          ++right_it;
        }
      }
      */
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

    VLMC_B_tree(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      Eigen::ArrayX4f cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      // Second pass - Comment this and select in dvstar.hpp row 120 to 128 to use old or new implementation.  
      container.second_pass(cached_context, background_order, offset_to_remove);
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
    ankerl::unordered_dense::map<int, RI_Kmer> container{};

  public: 
    VLMC_hashmap() = default;
    ~VLMC_hashmap() = default; 

    VLMC_hashmap(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      Eigen::ArrayX4f cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      // Second pass - Comment this and select in dvstar.hpp row 120 to 128 to use old or new implementation. 
      for (auto &[i_rep, kmer] : container) {
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        kmer.next_char_prob *= cached_context.row(offset).rsqrt();
      }
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
    int max_idx = 64; 
    std::array<RI_Kmer, 64> container_ibv{};
    std::vector<RI_Kmer> container_sorted{};  
    int max_kmer_index = 0;
    int min_kmer_index = 0;

  public: 
    VLMC_Combo() = default;
    ~VLMC_Combo() = default; 

    VLMC_Combo(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      Eigen::ArrayX4f cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      std::sort(container_sorted.begin(), container_sorted.end());

      // Second pass - Comment this and select in dvstar.hpp row 120 to 128 to use old or new implementation. 
      for (size_t i = 0; i < max_idx; i++){
        RI_Kmer kmer = container_ibv[i];
        if(kmer.is_null) continue;
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        container_ibv[i].next_char_prob *= cached_context.row(offset).rsqrt();
      }
      for (size_t i = 0; i < size(); i++){
        auto kmer = get(i);
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        get(i).next_char_prob *= cached_context.row(offset).rsqrt();
      }
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

class VLMC_Veb : public VLMC_Container {

  private: 
    veb::Veb_tree veb = veb::Veb_tree(500000);
    int min_index = INT_MAX;
    int max_index = -1;

  public: 
    VLMC_Veb() = default;
    ~VLMC_Veb() = default; 

    VLMC_Veb(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      
      std::ifstream ifs(path_to_bintree, std::ios::binary);
      cereal::BinaryInputArchive archive(ifs);

      Kmer input_kmer{};

      Eigen::ArrayX4f cached_context((int)std::pow(4, background_order), 4);
      std::vector<RI_Kmer> tmp_container{};

      auto offset_to_remove = 0;
      for (int i = 0; i < background_order; i++){
        offset_to_remove += std::pow(4, i); 
      }

      while (ifs.peek() != EOF){
        archive(input_kmer);
        RI_Kmer ri_kmer{input_kmer};

        if(input_kmer.length <= background_order){
          if (input_kmer.length + 1 > background_order){
            int offset = ri_kmer.integer_rep - offset_to_remove; 
            cached_context.row(offset) = ri_kmer.next_char_prob;
          }
        } else {
          if(ri_kmer.integer_rep > max_index) max_index = ri_kmer.integer_rep;
          if(ri_kmer.integer_rep < min_index) min_index = ri_kmer.integer_rep;
          tmp_container.push_back(ri_kmer);
        }
      }

      ifs.close();
      for(auto kmer : tmp_container){
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        kmer.next_char_prob *= cached_context.row(offset).rsqrt();
        veb::insert(veb, kmer);
      }
    } 

    size_t size() const override { return veb.size; }

    void push(const RI_Kmer &kmer) override { 
      if(veb.size < kmer.integer_rep){
        std::cout << "Too enourmous kmer : " << kmer.integer_rep << std::endl;
        std::cout << "Supports max : " << veb.size << std::endl;
        return;
      }
      veb::insert(veb, kmer); 
    }

    RI_Kmer &get(const int i) override { return null_kmer; }

    int get_max_kmer_index() const override { return max_index; }
    int get_min_kmer_index() const override { return min_index; }

    RI_Kmer find(const int i_rep) override { return veb::find(veb, i_rep); }

    void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f) override {
      RI_Kmer left_kmer = left_kmers.find(this->min_index);
      RI_Kmer right_kmer = right_kmers.find(this->min_index);
      while(true){
        // Check if right_kmers has this succeeding left_kmer
        // if, apply f
        if(left_kmer == right_kmer){
          f(left_kmer, right_kmer);
        }
        // Iterating left
        left_kmer = veb::succ(veb, left_kmer);
        if(left_kmer.integer_rep == -1){
          return;
        }
        right_kmer = right_kmers.find(left_kmer.integer_rep);
      }
    }
};

class VLMC_Set : public VLMC_Container {

  private:
    std::set<RI_Kmer> container{};
    int min_kmer = INT_MAX;
    int max_kmer = -1;

  public: 
    VLMC_Set() = default;
    ~VLMC_Set() = default; 

    VLMC_Set(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      
      std::ifstream ifs(path_to_bintree, std::ios::binary);
      cereal::BinaryInputArchive archive(ifs);

      Kmer input_kmer{};

      Eigen::ArrayX4f cached_context((int)std::pow(4, background_order), 4);
      std::vector<RI_Kmer> tmp_container{};

      auto offset_to_remove = 0;
      for (int i = 0; i < background_order; i++){
        offset_to_remove += std::pow(4, i); 
      }

      while (ifs.peek() != EOF){
        archive(input_kmer);
        RI_Kmer ri_kmer{input_kmer};

        if(input_kmer.length <= background_order){
          if (input_kmer.length + 1 > background_order){
            int offset = ri_kmer.integer_rep - offset_to_remove; 
            cached_context.row(offset) = ri_kmer.next_char_prob;
          }
        } else {
          if(ri_kmer.integer_rep > max_kmer) max_kmer = ri_kmer.integer_rep;
          if(ri_kmer.integer_rep < min_kmer) min_kmer = ri_kmer.integer_rep;
          tmp_container.push_back(ri_kmer);
        }
      }

      ifs.close();
      for(auto kmer : tmp_container){
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        kmer.next_char_prob *= cached_context.row(offset).rsqrt();
        container.insert(kmer);
      }
    }

    size_t size() const override { return container.size(); }

    void push(const RI_Kmer &kmer) override { container.insert(kmer); }

    RI_Kmer &get(const int i) override { return null_kmer; }

    int get_max_kmer_index() const override { return max_kmer; }
    int get_min_kmer_index() const override { return min_kmer; }

    std::set<RI_Kmer>::iterator begin_set(){
      return container.begin();
    }

    std::set<RI_Kmer>::iterator end_set(){
      return container.end();
    }

    RI_Kmer find(const int i_rep) override { return null_kmer; }

    void iterate_kmers(VLMC_Container &left_kmers, VLMC_Container &right_kmers,
    const std::function<void(const RI_Kmer &left, const RI_Kmer &right)> &f) override {
      auto left_kmers_c = static_cast<VLMC_Set&>(left_kmers);
      auto right_kmers_c = static_cast<VLMC_Set&>(right_kmers);
      std::set<RI_Kmer>::iterator l_it = left_kmers_c.begin_set();
      std::set<RI_Kmer>::iterator r_it = right_kmers_c.begin_set();
      while(l_it != left_kmers_c.end_set() && r_it != right_kmers_c.end_set()){
        if(*l_it == *r_it){
          f(*l_it, *r_it);
          l_it++;
          r_it++;
        } else if (*l_it < *r_it) {
          l_it++;
        } else {
          r_it++;
        }
      }
    }

};
}