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

class Kmer_Cluster {

  private: 
    // std::map<int, std::vector<Kmer_Pair>> container{}; 
    std::unordered_map<int, std::vector<Kmer_Pair>> container{}; 
    // std::unordered_multimap<int, Kmer_Pair> container{};
    size_t vlmc_count = 0; 

  public: 
    Kmer_Cluster() = default;
    ~Kmer_Cluster() = default; 

    int size(){ return vlmc_count; }

    void set_size(size_t count){ this->vlmc_count = count; }

    void push(const Kmer_Pair kmer_pair) { 
      container[kmer_pair.kmer.integer_rep].push_back(kmer_pair); 
    }

    void push_all(Kmer_Cluster cluster){
      auto begin_it = cluster.get_begin();
      auto end_it = cluster.get_end();
      while(begin_it != end_it){
        container[begin_it->first].insert(container[begin_it->first].end(), begin_it->second.begin(), begin_it->second.end());
        begin_it++; 
      }
    }

    std::vector<Kmer_Pair> &get(int bucket_num){
      return container[bucket_num]; 
    }

    int bucket_count() { return container.bucket_count(); }

    int experimental_bucket_count() { 
      auto count = 0;
      for (auto it = container.begin(); it != container.end(); it++){
        count++; 
      }
      return count; 
    }

    std::unordered_map<int, std::vector<Kmer_Pair>>::iterator get_begin(){
      return container.begin();
    }

    std::unordered_map<int, std::vector<Kmer_Pair>>::iterator get_end(){
      return container.end();
    }

    std::unordered_map<int, std::vector<Kmer_Pair>>::iterator find(int bucket_num){
      return container.find(bucket_num); 
    }

    int get_bucket(const Kmer_Pair& kmer_pair){
      return container.bucket(kmer_pair.id);
    }

    bool is_bucket_empty(int bucket_num){
      if (container.find(bucket_num) != container.end()){
        return true;
      }
      return false; 
      // return container.contains(bucket_num); 
    }

    void prettyPrint(){
      for (auto it = container.begin(); it != container.end(); it++){
      // for (int bucket_num = 0; bucket_num < container.max_size(); bucket_num++){
        std::cout << "Bucket number " << it->first; 
        if (it->second.size() < 1){
          std::cout << "()" << std::endl; 
          continue; 
        }
        for (auto it : it->second){ // (auto it = container.begin(bucket_num); it != container.end(bucket_num); ++it){
          std::cout << "(" << it.id << ", " << it.kmer.integer_rep << ") ";  
        }
        std::cout << std::endl;
      }
    }
};
}