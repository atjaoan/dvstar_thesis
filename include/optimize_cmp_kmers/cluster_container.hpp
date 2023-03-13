#pragma once 

#include <functional>
#include <filesystem>
#include <map>
#include <unordered_map>
 
#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_container.hpp"


/*
  Stores multiple VLMCs in a container. 
*/
namespace container {

template <typename VC> 
class Cluster_Container {

  private: 
    std::vector<VC> container{};

  public: 
    Cluster_Container() = default;
    ~Cluster_Container() = default; 

    Cluster_Container(const size_t i) : container(i) { }

    size_t size() const { return container.size(); }

    void push(const VC vlmc) { container.push_back(vlmc); }

    void for_each(const std::function<void(VC vlmc)> &f) {
      for (auto vlmc : container){
        f(vlmc);
      }
    }

    VLMC_Container &get(const int i) { return container[i]; }

    VC &operator[](size_t index){ return container[index]; }

    const VC &operator[](size_t index) const { return container[index]; }

};

struct Kmer_Pair {
  RI_Kmer kmer;
  size_t id; 

  Kmer_Pair(RI_Kmer kmer, size_t id){
    this->kmer = kmer; 
    this->id = id; 
  }

  Kmer_Pair() = default;

  ~Kmer_Pair() = default; 
};

using bucket_t = std::unordered_multimap<int, Kmer_Pair>::local_iterator;

class Kmer_Cluster {

  private: 
    // std::map<int, std::vector<Kmer_Pair>> container{}; 
    std::unordered_multimap<int, Kmer_Pair> container{};
    size_t vlmc_count = 0; 

  public: 
    Kmer_Cluster() = default;
    ~Kmer_Cluster() = default; 

    int size(){ return vlmc_count; }

    void set_size(size_t count){ this->vlmc_count = count; }

    void push(const Kmer_Pair kmer_pair) { 
      container.insert({kmer_pair.kmer.integer_rep, kmer_pair}); 
    }

    std::vector<Kmer_Pair> get(int bucket_num){
      std::vector<Kmer_Pair> returning_vector{};
      for( auto it = container.begin(bucket_num); it!= container.end(bucket_num); ++it){
        returning_vector.push_back(it->second);
      } 
      return returning_vector;
    }

    bucket_t get_bucket_begin(int bucket_num){
      return container.begin(bucket_num);
    }

    bucket_t get_bucket_end(int bucket_num){
      return container.end(bucket_num);
    }

    int bucket_count() { return container.bucket_count(); }

    int get_bucket(const Kmer_Pair& kmer_pair){
      return container.bucket(kmer_pair.id);
    }

    bool is_bucket_empty(int bucket_num){
      return container.bucket_size(bucket_num) == 0; 
    }

    void prettyPrint(){
      size_t print_size = 0; 
      for (int bucket_num = 0; bucket_num < container.bucket_count(); bucket_num++){
        if (print_size > 10){
          break; 
        }
        if (container.bucket_size(bucket_num) < 2){
          continue; 
        }
        for (auto it = container.begin(bucket_num); it != container.end(bucket_num); ++it){
          std::cout << "(" << it->first << ", {" << it->second.id << ", " << it->second.kmer.next_char_prob[0] << "}) ";  
        }
        std::cout << std::endl;
        print_size++; 
      }
    }
};
}