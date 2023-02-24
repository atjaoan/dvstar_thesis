#pragma once

#include <filesystem>
#include <stdlib.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

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
      cluster[index] = VC{paths[index]};
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
}
