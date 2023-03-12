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

#include "cluster_container.hpp"
#include "vlmc_container.hpp"
#include "parallel.hpp"

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;
using RI_Kmer = container::RI_Kmer;

namespace cluster{
  
template <typename VC>  
container::Cluster_Container<VC> get_cluster(const std::filesystem::path &directory, const size_t nr_cores_to_use, const size_t background_order){
  std::vector<std::filesystem::path> paths{};

  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    paths.push_back(dir_entry.path());
  }

  container::Cluster_Container<VC> cluster{paths.size()};

  auto fun = [&](size_t start_index, size_t stop_index) {
    for (int index = start_index; index < stop_index; index++){
      cluster[index] = VC(paths[index], background_order);
    }
  };

  parallel::parallelize(paths.size(), fun, nr_cores_to_use);

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

container::Kmer_Cluster get_kmer_cluster(const std::filesystem::path &directory, const size_t background_order = 0){
  container::Kmer_Cluster cluster{};
  size_t id = 0; 
  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    std::ifstream ifs(dir_entry.path(), std::ios::binary);
    cereal::BinaryInputArchive archive(ifs);
    Kmer input_kmer{};
    std::vector<container::RI_Kmer> input_vector{}; 

    Eigen::ArrayX4d cached_context((int)std::pow(4, background_order), 4); 

    auto offset_to_remove = 0;
    for (int i = 0; i < background_order; i++){
      offset_to_remove += std::pow(4, i); 
    }

    while (ifs.peek() != EOF){
      archive(input_kmer);
      RI_Kmer ri_kmer{input_kmer};
      if(ri_kmer.length <= background_order){
        if (ri_kmer.length + 1 > background_order){
          int offset = ri_kmer.integer_rep - offset_to_remove; 
          cached_context.row(offset) = ri_kmer.next_char_prob;
        }
      } else {
        input_vector.push_back(ri_kmer); 
      }
    }
    ifs.close();

    for (RI_Kmer kmer : input_vector){
      int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
      int offset = background_idx - offset_to_remove;
      kmer.next_char_prob *= cached_context.row(offset).rsqrt();
      cluster.push(container::Kmer_Pair{kmer, id}); 
    }
    id++; 
  }

  cluster.set_size(id); 

  return cluster; 
}
}
