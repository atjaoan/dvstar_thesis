#pragma once

#include <filesystem>
#include <stdlib.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>
#include <functional>
#include <limits.h>
#include <exception>
#include <algorithm>
#include <mutex>

#include "cluster_container.hpp"
#include "vlmc_container.hpp"
#include "parallel.hpp"
#include "utils.hpp"

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;
using RI_Kmer = container::RI_Kmer;

namespace cluster{
  
template <typename VC>  
container::Cluster_Container<VC> get_cluster(const std::filesystem::path &directory, size_t nr_cores_to_use, 
      const size_t background_order, const int set_size = -1){
  std::vector<std::filesystem::path> paths{};

  if(nr_cores_to_use > 4) nr_cores_to_use = 4;
  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    paths.push_back(dir_entry.path());
  }

  size_t paths_size = paths.size();
  if ((set_size != -1) && (set_size < paths_size)){
    paths_size = set_size; 
  }

  container::Cluster_Container<VC> cluster{paths_size};

  auto fun = [&](size_t start_index, size_t stop_index) {
    for (int index = start_index; index < stop_index; index++){
      cluster[index] = VC(paths[index], background_order);
    }
  };

  parallel::parallelize(paths_size, fun, nr_cores_to_use);

  return cluster; 
}

/*
  Will keep this non-parallel implementation here for testing purposes
*/
template <typename VC>  
container::Cluster_Container<VC> old_get_cluster(const std::filesystem::path &directory){
  container::Cluster_Container<VC> cluster{};
  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    cluster.push(VC{dir_entry.path()});
  }
  return cluster; 
}

// Return as reference?
container::Kmer_Cluster get_kmer_cluster(const std::filesystem::path &directory, const size_t background_order = 0){
  std::vector<std::filesystem::path> paths{};

  auto nr_cores_to_use = 4; 
  auto set_size = -1; 

  if(nr_cores_to_use > 4) nr_cores_to_use = 4;
  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    paths.push_back(dir_entry.path());
  }

  size_t paths_size = paths.size();
  if ((set_size != -1) && (set_size < paths_size)){
    paths_size = set_size; 
  }

  container::Kmer_Cluster final_cluster{}; 

  std::mutex cluster_lock;

  auto fun = [&](size_t start_index, size_t stop_index) {
    container::Kmer_Cluster cluster{}; 
    for (int index = start_index; index < stop_index; index++){
      std::vector<container::RI_Kmer> input_vector{};  
      Eigen::ArrayX4f cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { input_vector.push_back(kmer); }; 

      int offset_to_remove = container::load_VLMCs_from_file(paths[index], cached_context, fun, background_order); 

      for (RI_Kmer kmer : input_vector){
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        kmer.next_char_prob *= cached_context.row(offset).rsqrt();
        cluster.push(container::Kmer_Pair{kmer, index}); 
      }
    } 

    std::lock_guard<std::mutex> lock(cluster_lock);
    final_cluster.push_all(cluster); 
  };

  parallel::parallelize(paths_size, fun, nr_cores_to_use);

  final_cluster.set_size(paths_size); 


  return final_cluster; 
}
}
