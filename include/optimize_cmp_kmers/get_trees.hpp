#pragma once

#include <filesystem>

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;

namespace get_trees{
  
template<typename CC, typename VC>
void get_trees(const std::filesystem::path &directory, CC trees){

  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    VC tree(dir_entry.path());
    trees.push(tree);
  }

  std::cout << trees.size() << std::endl; 
}
}
