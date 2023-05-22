#pragma once

#include <filesystem>
#include <thread>
#include <mutex>

#include "cluster_container.hpp"
#include "global_aliases.hpp"
#include "parallel.hpp"

namespace cluster{

using RI_Kmer = container::RI_Kmer;
  
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

std::vector<container::Kmer_Cluster> get_kmer_cluster(const std::filesystem::path &directory, BS::thread_pool& pool, 
          const size_t background_order = 0, const int set_size = -1){
  std::vector<std::filesystem::path> paths{};

  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    paths.push_back(dir_entry.path());
  }

  size_t paths_size = paths.size();
  if ((set_size != -1) && (set_size < paths_size)){
    paths_size = set_size; 
  }

  std::vector<container::Kmer_Cluster> clusters(std::ceil(paths_size / std::floor(std::sqrt(paths_size)))); 

  auto fun = [&](size_t start_index, size_t stop_index, size_t idx) {
    for (int index = start_index; index < stop_index; index++){
      std::vector<container::RI_Kmer> input_vector{};  
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer &kmer) { input_vector.push_back(kmer); }; 

      int offset_to_remove = container::load_VLMCs_from_file(paths[index], cached_context, fun, background_order); 

      for (RI_Kmer kmer : input_vector){
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        kmer.next_char_prob *= cached_context.row(offset).rsqrt();
        clusters[idx].push(container::Kmer_Pair{kmer, index - start_index}); 
      }
    } 
    clusters[idx].set_size(stop_index - start_index);
  };

  parallel::parallelize_kmer_major(paths_size, fun, pool);

  return clusters; 
}
}
