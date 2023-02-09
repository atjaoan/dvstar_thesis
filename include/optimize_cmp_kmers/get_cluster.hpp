#pragma once

#include <filesystem>

#include "cluster_container.hpp"
#include "vlmc_container.hpp"

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;

namespace cluster{
  
template <typename VC>  
void get_cluster(const std::filesystem::path &directory, container::Cluster_Container &cluster){
  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    VC vlmc{dir_entry.path()};
    cluster.push(std::make_shared<VC>(vlmc));
  }
}
}
