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

container::Kmer_Cluster get_kmer_cluster(const std::filesystem::path &directory){
  container::Kmer_Cluster cluster{};
  size_t id = 0; 
  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    std::ifstream ifs(dir_entry.path(), std::ios::binary);
    cereal::BinaryInputArchive archive(ifs);
    Kmer input_kmer{};
    
    while (ifs.peek() != EOF){
      archive(input_kmer);
      container::RI_Kmer ri_kmer{input_kmer};
      container::Kmer_Pair kmer_pair{ri_kmer, id};
      cluster.push(kmer_pair);
    }
    ifs.close();
    id++; 
  }
  return cluster; 
}



}
