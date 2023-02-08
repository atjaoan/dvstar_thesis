#pragma once

#include <filesystem>

#include "cluster_container.hpp"
#include "vlmc_container.hpp"

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;

namespace get_trees{
  
template <typename VC>  
void get_trees(const std::filesystem::path &directory, container::Cluster_Container trees){

  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    VC tree(dir_entry.path());
    std::cout << tree.size() << std::endl; 
    trees.push(tree);
  }

  std::cout << trees.size() << std::endl; 
}
}
