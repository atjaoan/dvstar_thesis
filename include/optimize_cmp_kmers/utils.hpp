#pragma once

#include <filesystem>
#include <thread>

#include "cluster_container.hpp"

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;
using RI_Kmer = container::RI_Kmer;
using matrix_t = Eigen::MatrixXd;

namespace utils {
  
void matrix_recursion(size_t start_index_left, size_t stop_index_left, size_t start_index_right, size_t stop_index_right, 
  const std::function<void(size_t &left, size_t &right)> &fun){
  auto diff_left = stop_index_left - start_index_left;
  auto diff_right = stop_index_right - start_index_right;
  if (diff_left == 1 && diff_right == 1){
    fun(start_index_left, start_index_right); 
  } else if (diff_right > diff_left){
    auto new_right_index = (stop_index_right + start_index_right) / 2;
    matrix_recursion(start_index_left, stop_index_left, start_index_right, new_right_index, fun);
    matrix_recursion(start_index_left, stop_index_left, new_right_index, stop_index_right, fun);
  } else {
    auto new_left_index = (stop_index_left + start_index_left) / 2;
    matrix_recursion(start_index_left, new_left_index, start_index_right, stop_index_right, fun);
    matrix_recursion(new_left_index, stop_index_left, start_index_right, stop_index_right, fun);
  }
}

void half_matrix_recursion(size_t start_index_left, size_t stop_index_left, size_t start_index_right, size_t stop_index_right, 
  const std::function<void(size_t &left, size_t &right)> &fun){
  if (start_index_left < start_index_right || stop_index_left <= stop_index_right){
    auto diff_left = stop_index_left - start_index_left;
    auto diff_right = stop_index_right - start_index_right; 
    if (diff_left == 1 && diff_right == 1){
      if (start_index_left != start_index_right){
        fun(start_index_left, start_index_right); 
      }
    } else if (diff_right > diff_left){
      auto new_right_index = (stop_index_right + start_index_right) / 2;
      half_matrix_recursion(start_index_left, stop_index_left, start_index_right, new_right_index, fun);
      half_matrix_recursion(start_index_left, stop_index_left, new_right_index, stop_index_right, fun);
    } else {
      auto new_left_index = (stop_index_left + start_index_left) / 2;
      half_matrix_recursion(start_index_left, new_left_index, start_index_right, stop_index_right, fun);
      half_matrix_recursion(new_left_index, stop_index_left, start_index_right, stop_index_right, fun);
    }
  } 
}

int load_VLMCs_from_file(const std::filesystem::path &path_to_bintree, Eigen::ArrayX4d &cached_context, 
        const std::function<void(const RI_Kmer &kmer)> f, const size_t background_order = 0) {
  std::ifstream ifs(path_to_bintree, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);
  Kmer input_kmer{};

  int offset_to_remove = 0;
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
      f(ri_kmer); 
    }
  }
  ifs.close();

  return offset_to_remove; 
}

int get_used_cores(size_t requested_cores, size_t size){
  const size_t processor_count = std::thread::hardware_concurrency();
  size_t used_cores = 1;
  if(requested_cores > size){
    used_cores = size;
  } else if(requested_cores <= processor_count){
      used_cores = requested_cores;
  } else {
    used_cores = processor_count;
  }
  return used_cores;
}

void print_matrix(matrix_t distance_matrix){
  for (size_t i = 0; i < distance_matrix.rows(); i++)
  {
    for (size_t j = 0; j < distance_matrix.cols(); j++)
    {
      std::cout << distance_matrix(i,j) << " ";
    }
    std::cout << std::endl;
  }
}

std::string get_filename(std::filesystem::path path){
  if (path.has_filename()){
    return path.parent_path().filename().u8string() + "_" + path.filename().u8string();
  } else {
    return path.parent_path().parent_path().filename().u8string() + "_" + path.parent_path().filename().u8string();
  }
}
}