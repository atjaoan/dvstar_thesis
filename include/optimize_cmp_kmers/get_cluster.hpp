#pragma once

#include <filesystem>

#include "cluster_container.hpp"
#include "vlmc_container.hpp"

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;

namespace cluster{
  
template <typename VC>  
container::Cluster_Container<VC> get_cluster(const std::filesystem::path &directory){
  container::Cluster_Container<VC> cluster{};
  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    cluster.push(VC{dir_entry.path()});
  }
  return cluster; 
}
}
