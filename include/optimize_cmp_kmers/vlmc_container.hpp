#pragma once 

#include <functional>
#include <filesystem>
#include <limits.h>
#include <exception>
#include <algorithm>
#include <set>
#include <execution>
#include <math.h>
#include <Eigen/Core> 

#include "vlmc_from_kmers/kmer.hpp"
#include "optimize_cmp_kmers/read_in_kmer.hpp"
#include "global_aliases.hpp"
#include "b_tree.hpp"
#include "robin_hood.h"
#include "unordered_dense.h"
//#include "veb_tree.hpp"
#include "veb_array.hpp"
#include "eytzinger_array.hpp"
#include "b_tree_alt.hpp"

/*
  Stores VLMC (multiple k-mers) in a container. 
*/
constexpr int misses_before_skip = 6;

using Kmer = vlmc::VLMCKmer; 

namespace container{

int load_VLMCs_from_file(const std::filesystem::path &path_to_bintree, eigenx_t &cached_context, 
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

out_t normalise_dvstar(out_t dot_product, out_t left_norm,
                        out_t right_norm) {

  left_norm = std::sqrt(left_norm);
  right_norm = std::sqrt(right_norm);
  if (left_norm == 0 || right_norm == 0) {
    return 1.0;
  } else {
    out_t Dvstar = dot_product / (left_norm * right_norm);

    out_t dvstar = 0.5 * (1 - Dvstar);

    out_t angular_distance = 2 * std::acos(Dvstar) / M_PI;
    if (isnan(angular_distance)) {
      return 0.0;
    } else {
      return angular_distance;
    }
  }
}

/*
  Storing Kmers in a unsorted vector.
*/
class VLMC_vector {

  private: 
    std::vector<RI_Kmer> container{}; 
    RI_Kmer null_kmer{};

  public: 
    VLMC_vector() = default;
    ~VLMC_vector() = default; 

    VLMC_vector(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      for (size_t i = 0; i <= get_max_kmer_index(); i++){
        RI_Kmer kmer = get(i);
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        get(i).next_char_prob *= cached_context.row(offset).rsqrt();
      }
    } 

    RI_Kmer* begin(){
      return &container[0];
    }

    RI_Kmer* end(){
      return begin() + container.size();
    }

    size_t size() const { return container.size(); }

    void push(const RI_Kmer &kmer) { container.push_back(kmer); }

    RI_Kmer &get(const int i) { return container[i]; }

    int get_max_kmer_index() const { return container.size() - 1; }
    int get_min_kmer_index() const { return 0; }

    RI_Kmer find(const int i_rep) {
      for (size_t i = 0; i < container.size(); i++){
        if (container[i].integer_rep==i_rep) {
          return container[i]; 
        }
      }
      return null_kmer; 
    }
};

out_t iterate_kmers(VLMC_vector &left_kmers, VLMC_vector &right_kmers) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;

  for (size_t i = left_kmers.get_min_kmer_index() ; i <= left_kmers.get_max_kmer_index(); i++) {
    const RI_Kmer &left_kmer = left_kmers.get(i);
    if (left_kmer.is_null){
      continue; 
    }
    auto right_kmer = right_kmers.find(left_kmer.integer_rep);
    if (!right_kmer.is_null){
      dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
      left_norm += left_kmer.next_char_prob.square().sum();
      right_norm += right_kmer.next_char_prob.square().sum();
    } 
  }
  return normalise_dvstar(dot_product, left_norm, right_norm);
}

/*
  Storing Kmers in a sorted vector.
*/
class VLMC_sorted_vector {

  private: 
    RI_Kmer null_kmer{};

  public: 
    std::vector<RI_Kmer> container{};
    VLMC_sorted_vector() = default;
    ~VLMC_sorted_vector() = default; 

    VLMC_sorted_vector(const std::filesystem::path &path_to_bintree, const size_t background_order = 0, bool use_new = false) {
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);
      
      std::sort(std::execution::seq, container.begin(), container.end());
      for (size_t i = 0; i < size(); i++){
        RI_Kmer kmer = get(i);
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        get(i).next_char_prob *= cached_context.row(offset).rsqrt();
      }
      //int k = 0;
      //for(auto kmer : tmp_container){
      //  std::cout << kmer.next_char_prob
      //  k++;
      //  if(k > 50) break;
      //}
    } 

    size_t size() const { return container.size(); }

    void push(const RI_Kmer &kmer) { container.push_back(kmer); }

    std::vector<RI_Kmer>::iterator begin() { return container.begin(); };
    std::vector<RI_Kmer>::iterator end() { return container.end(); };

    RI_Kmer &get(const int i) { return container[i]; }

    int get_max_kmer_index() const { return container.size() - 1; }
    int get_min_kmer_index() const { return 0; }

    RI_Kmer find(const int i_rep) {
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
};

out_t iterate_kmers(VLMC_sorted_vector &left_kmers, VLMC_sorted_vector &right_kmers) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;
  //int f_applied = 0;
  auto right_it = right_kmers.begin();
  auto right_end = right_kmers.end();
  auto left_it = left_kmers.begin();
  auto left_end = left_kmers.end();
  
  while(left_it != left_end && right_it != right_end){
    auto left_kmer = *left_it;
    auto right_kmer = *right_it;
    if(left_kmer == right_kmer){
      dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
      left_norm += left_kmer.next_char_prob.square().sum();
      right_norm += right_kmer.next_char_prob.square().sum();
      ++left_it;
      ++right_it;
      //f_applied++;
    } else if(left_kmer < right_kmer) {
      ++left_it;
    }
    else ++right_it;
  }
  //std::cout << f_applied << "\n";
  return normalise_dvstar(dot_product, left_norm, right_norm); 
}

/*
  Storing Kmers in a B-tree.
*/

class VLMC_B_tree {

  private: 
    RI_Kmer null_kmer{};

  public: 
    b_tree::BTree container{3};
    VLMC_B_tree() = default;
    ~VLMC_B_tree() = default; 

    VLMC_B_tree(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);
 
      container.second_pass(cached_context, background_order, offset_to_remove);
    } 

    size_t size() const { return 0; }

    void push(const RI_Kmer &kmer) { container.insert(kmer); }

    RI_Kmer &get(const int i) { return null_kmer; }

    int get_max_kmer_index() const { return -1; }
    int get_min_kmer_index() const { return 0; }

    RI_Kmer find(const int i_rep) {
      return container.search(i_rep); 
    }
};

out_t iterate_kmers(VLMC_B_tree &left_kmers, VLMC_B_tree &right_kmers) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;

  left_kmers.container.for_each([&](const RI_Kmer &left_kmer) {
    RI_Kmer right_kmer = right_kmers.find(left_kmer.integer_rep);
    if (!right_kmer.is_null){
      dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
      left_norm += left_kmer.next_char_prob.square().sum();
      right_norm += right_kmer.next_char_prob.square().sum();
    } 
  });

  return normalise_dvstar(dot_product, left_norm, right_norm);
}

/*
  Storing Kmers in a unordered map (HashMap).
*/
class VLMC_hashmap {

  private: 
    RI_Kmer null_kmer{};

  public: 
    ankerl::unordered_dense::map<int, RI_Kmer> container{};
    VLMC_hashmap() = default;
    ~VLMC_hashmap() = default; 

    VLMC_hashmap(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      for (auto &[i_rep, kmer] : container) {
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        kmer.next_char_prob *= cached_context.row(offset).rsqrt();
      }
    } 

    size_t size() const { return container.size(); }

    void push(const RI_Kmer &kmer) { container[kmer.integer_rep] = kmer; }

    RI_Kmer &get(const int i) { return container[i]; }

    int get_max_kmer_index() const { return container.size() - 1; }
    int get_min_kmer_index() const { return 0; }

    RI_Kmer find(const int i_rep) {
      auto res = container.find(i_rep);
      if (res != container.end()){
        return res->second; 
      }
      return null_kmer; 
    }
};

out_t iterate_kmers(VLMC_hashmap &left_kmers, VLMC_hashmap &right_kmers) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;

  for (auto &[i_rep, left_kmer] : left_kmers.container) {
    auto right_kmer = right_kmers.find(i_rep); 
    if (!right_kmer.is_null) {
      dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
      left_norm += left_kmer.next_char_prob.square().sum();
      right_norm += right_kmer.next_char_prob.square().sum();
    } 
  }

  return normalise_dvstar(dot_product, left_norm, right_norm);
}

class VLMC_Veb {

  private: 
    int min_index = INT_MAX;
    int max_index = -1;
    RI_Kmer null_kmer{};

  public:
    veb::Veb_array *veb;
    VLMC_Veb() = default;
    ~VLMC_Veb() = default; 

    VLMC_Veb(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto tmp_container = std::vector<RI_Kmer>{};
      auto fun = [&](const RI_Kmer &kmer) { tmp_container.push_back(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);
      
      std::sort(std::execution::seq, tmp_container.begin(), tmp_container.end());
      for (auto &kmer : tmp_container){
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        kmer.next_char_prob *= cached_context.row(offset).rsqrt();
      }
      veb = new veb::Veb_array(tmp_container);
    }

    size_t size() const { return veb->n + 1; }

    RI_Kmer &get(const int i) {
      return veb->get_from_array(i);
    }

    int get_max_kmer_index() const { return max_index; }
    int get_min_kmer_index() const { return min_index; }

    RI_Kmer find(const int i_rep) { return null_kmer; }
};

out_t iterate_kmers(VLMC_Veb &left_kmers, VLMC_Veb &right_kmers) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;
  int i = 0;
  while(i < left_kmers.veb->n){
    RI_Kmer& left_kmer = left_kmers.veb->a[i];
    RI_Kmer& right_kmer = right_kmers.get(left_kmer.integer_rep);
    if(left_kmer == right_kmer){
      dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
      left_norm += left_kmer.next_char_prob.square().sum();
      right_norm += right_kmer.next_char_prob.square().sum();
    }
    i++;
  }
  return normalise_dvstar(dot_product, left_norm, right_norm);
}

class VLMC_Eytzinger {

  private: 
    int min_index = INT_MAX;
    int max_index = -1;

  public: 
    array::Ey_array *arr;
    VLMC_Eytzinger() = default;
    ~VLMC_Eytzinger() = default; 

    VLMC_Eytzinger(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto tmp_container = std::vector<RI_Kmer>{};
      auto fun = [&](const RI_Kmer &kmer) { tmp_container.push_back(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);
      
      std::sort(std::execution::seq, tmp_container.begin(), tmp_container.end());
      for (auto &kmer : tmp_container){
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        kmer.next_char_prob *= cached_context.row(offset).rsqrt();
      }
      arr = new array::Ey_array(tmp_container);
    } 

    size_t size() const { return arr->size + 1; }

    RI_Kmer &get(const int i) { ;
      return arr->get_from_array(i);
    }

    int get_max_kmer_index() const { return max_index; }
    int get_min_kmer_index() const { return min_index; }
};

out_t iterate_kmers(VLMC_Eytzinger &left_kmers, VLMC_Eytzinger &right_kmers) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;
  int i = 0;
  while(i <= left_kmers.arr->size){
    RI_Kmer& left_kmer = left_kmers.arr->ey_sorted_kmers[i];
    RI_Kmer& right_kmer = right_kmers.get(left_kmer.integer_rep);
    if(left_kmer == right_kmer){
      dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
      left_norm += left_kmer.next_char_prob.square().sum();
      right_norm += right_kmer.next_char_prob.square().sum();
    }
    i++;
  }
  return normalise_dvstar(dot_product, left_norm, right_norm);
}

class VLMC_Alt_Btree {

  private: 
    int min_index = INT_MAX;
    int max_index = -1;

  public: 
    array::B_Tree *arr;
    VLMC_Alt_Btree() = default;
    ~VLMC_Alt_Btree() = default; 

    VLMC_Alt_Btree(const std::filesystem::path &path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto tmp_container = std::vector<RI_Kmer>{};
      auto fun = [&](const RI_Kmer &kmer) { tmp_container.push_back(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);
      
      std::sort(std::execution::seq, tmp_container.begin(), tmp_container.end());
      for (auto &kmer : tmp_container){
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        kmer.next_char_prob *= cached_context.row(offset).rsqrt();
      }
      arr = new array::B_Tree(tmp_container);
    } 

    size_t size() const { return arr->size + 1; }

    RI_Kmer &get(const int i) { ;
      return arr->get_from_array(i);
    }

    int get_max_kmer_index() const { return max_index; }
    int get_min_kmer_index() const { return min_index; }
}; 

out_t iterate_kmers(VLMC_Alt_Btree &left_kmers, VLMC_Alt_Btree &right_kmers) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;
  int i = 0;
  while(i < left_kmers.arr->size){
    RI_Kmer& left_kmer = left_kmers.arr->a[i];
    RI_Kmer& right_kmer = right_kmers.get(left_kmer.integer_rep);
    if(left_kmer == right_kmer){
      dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
      left_norm += left_kmer.next_char_prob.square().sum();
      right_norm += right_kmer.next_char_prob.square().sum();
    }
    i++;
  }
  return normalise_dvstar(dot_product, left_norm, right_norm);
}

/*
  Storing Kmers in a sorted vector with a summary structure to skip past misses.
*/

struct Min_max_node {
  RI_Kmer* block_start;
  int min, max; 

  Min_max_node(RI_Kmer* kmer, int min, int max){
    this->block_start = kmer; 
    this->min = min;
    this->max = max;
  }

  Min_max_node() = default;
  ~Min_max_node() = default; 
};

class VLMC_sorted_search {

  private: 
    RI_Kmer null_kmer{};
    int skip_size;

  public: 
    std::vector<RI_Kmer> container{};
    std::vector<Min_max_node> summary{};
    int place_in_summary = 0;
    VLMC_sorted_search() = default;
    ~VLMC_sorted_search() = default; 

    VLMC_sorted_search(const std::filesystem::path &path_to_bintree, const size_t background_order = 0, bool use_new = false) {
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { push(kmer); }; 

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);
      
      std::sort(std::execution::seq, container.begin(), container.end());
      for (size_t i = 0; i < size(); i++){
        RI_Kmer kmer = get(i);
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        get(i).next_char_prob *= cached_context.row(offset).rsqrt();
      }
      // Build summary
      skip_size = std::sqrt(container.size());
      summary.reserve(skip_size);
      //std::cout << "size: " << container.size() << " skip size: " << skip_size << "\n";
      int i = 0;
      for(;i < container.size() - skip_size; i += skip_size){
        summary.push_back(Min_max_node(&container[i], container[i].integer_rep, container[i+skip_size - 1].integer_rep));
      }
      summary.push_back(Min_max_node(&container[i], container[i].integer_rep, container[size() - 1].integer_rep));
      for(auto& node : summary){
        //std::cout << node.block_start << " - " << node.min << " - " << node.max << "\n";
      }
    } 

    size_t size() const { return container.size(); }

    void push(const RI_Kmer &kmer) { container.push_back(kmer); }

    std::vector<RI_Kmer>::iterator begin() { return container.begin(); };
    std::vector<RI_Kmer>::iterator end() { return container.end(); };

    RI_Kmer &get(const int i) { return container[i]; }

    int get_max_kmer_index() const { return container.size() - 1; }
    int get_min_kmer_index() const { return 0; }

    RI_Kmer find(const int i_rep) {
      return null_kmer;
    }

    RI_Kmer* find_block_start(int i_rep){
      for(int i = place_in_summary; i < skip_size; i++){
        if(i_rep >= summary[i].min && i_rep < summary[i].max){
          place_in_summary = i;
          return summary[i].block_start;
        }
      }
      return nullptr;
    }
};

out_t iterate_kmers(VLMC_sorted_search &left_kmers, VLMC_sorted_search &right_kmers) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;

  //int f_applied = 0;

  RI_Kmer* left_kmer = &left_kmers.container[0];
  RI_Kmer* right_kmer = &right_kmers.container[0];
  left_kmers.place_in_summary = 0;
  right_kmers.place_in_summary = 0;
  auto left_end = &left_kmers.container[left_kmers.container.size()];
  auto right_end = &right_kmers.container[right_kmers.container.size()];
  
  while(left_kmer < left_end && right_kmer < right_end){
    if(*left_kmer == *right_kmer){
      dot_product += (left_kmer->next_char_prob * right_kmer->next_char_prob).sum();
      left_norm += left_kmer->next_char_prob.square().sum();
      right_norm += right_kmer->next_char_prob.square().sum();
      ++left_kmer;
      ++right_kmer;
      //f_applied++;
    } else if(*left_kmer < *right_kmer) {
      if(left_kmers.summary[left_kmers.place_in_summary].max < right_kmer->integer_rep){
        left_kmer = left_kmers.find_block_start(right_kmer->integer_rep);
      } else {
        ++left_kmer;
      }
    }
    else {
      if(right_kmers.summary[right_kmers.place_in_summary].max < left_kmer->integer_rep){
        right_kmer = right_kmers.find_block_start(left_kmer->integer_rep);
      } else {
        ++right_kmer;
      }
    }
    if(left_kmer == nullptr || right_kmer == nullptr) break;
  }
  //std::cout << f_applied << "\n";
  return normalise_dvstar(dot_product, left_norm, right_norm); 
}

}