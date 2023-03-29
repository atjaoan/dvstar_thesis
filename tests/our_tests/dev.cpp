#pragma once 

#include <functional>
#include <filesystem>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
#include <math.h>
#include <iostream>
#include <string>

#include "read_in_kmer.hpp"

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;
using RI_Kmer = container::RI_Kmer; 
using matrix_t = Eigen::MatrixXf;
#include "utils.hpp"

float normalise_dvstar_test(float dot_product, float left_norm,
                        float right_norm) {

  left_norm = std::sqrt(left_norm);
  right_norm = std::sqrt(right_norm);
  if (left_norm == 0 || right_norm == 0) {
    return 1.0;
  } else {
    float Dvstar = dot_product / (left_norm * right_norm);

    float dvstar = 0.5 * (1 - Dvstar);

    float angular_distance = 2 * std::acos(Dvstar) / M_PI;
    if (isnan(angular_distance)) {
      return 0.0;
    } else {
      return angular_distance;
    }
  }
}

void get_kmer_vector_dev(std::filesystem::path path, std::vector<std::vector<RI_Kmer>>& vlmcs, int max_size, int offset){
  std::ifstream ifs(path, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);
  std::vector<RI_Kmer> vlmc{};
  vlmc::VLMCKmer input_kmer{};
  int curr_size = 0;
  while (ifs.peek() != EOF && curr_size < max_size){
    archive(input_kmer);
    RI_Kmer ri_kmer(input_kmer); 
    vlmc.push_back(ri_kmer);
    curr_size++;
  }
  //curr_size = 0;
  //ifs.seekg(0, ifs.beg);
  ifs.close();
  std::sort(vlmc.begin(), vlmc.end());
  vlmc.reserve(vlmc.size() + offset);
  vlmcs.push_back(vlmc);
}

matrix_t iterate_kmers_bench_dev(int max_size){
  std::filesystem::path path_fst{"./data/human_VLMCs"};
  std::vector<std::vector<RI_Kmer>> vlmcs{};
  int offset = 0;
  for (const auto& dir_entry : recursive_directory_iterator(path_fst)) {
    get_kmer_vector_dev(dir_entry.path(), vlmcs, max_size, offset);
    offset += 257;
  }
  
  matrix_t distances{vlmcs.size(), vlmcs.size()};
  auto start = std::chrono::steady_clock::now();

  auto fun = [&](auto x, auto y) {
    std::vector<RI_Kmer> &left_vlmc = vlmcs[x];
    std::vector<RI_Kmer> &right_vlmc = vlmcs[y];

    float dot_product = 0.0;
    float left_norm = 0.0;
    float right_norm = 0.0; 

    auto dvstar_fun = [&](auto &left_v, auto &right_v) {
      dot_product += (left_v.next_char_prob * right_v.next_char_prob).sum();
      left_norm += left_v.next_char_prob.square().sum();
      right_norm += right_v.next_char_prob.square().sum();
    };

    auto right_it = right_vlmc.begin();
    auto right_end = right_vlmc.end();
    for(auto &left_kmer : left_vlmc){
      while((*right_it) < left_kmer){
        ++right_it;
        if(right_it == right_end) break;
      }
      if(left_kmer == (*right_it)){
        dvstar_fun(left_kmer, *right_it);
        ++right_it;
      }
    }

    distances(x,y) = normalise_dvstar_test(dot_product, left_norm, right_norm); 
  };

  utils::matrix_recursion(0, vlmcs.size(), 0, vlmcs.size(), fun); 

  auto end = std::chrono::steady_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  std::cout << "Time : " << time << " [micro sec]" << std::endl; 
  for(int i = 0; i < vlmcs.size(); i++){
    std::cout << "Size of vectors: " << sizeof(vlmcs[i]) + sizeof(RI_Kmer) * vlmcs[i].capacity() << " bytes" << "\n";
  }
  return distances;
}

matrix_t iterate_vectors_int(int max_size){
  auto container = std::vector<int>(max_size);
  matrix_t distances{24, 24};
  
  auto container_it = container.begin();
  auto container_end = container.end();
  int i = 0;
  while(container_it != container_end){
    //*container_it += 2;
    distances(i%24,i%24) = (*container_it) + 5;
    ++container_it;
    i++;
  }

  return distances;
}

matrix_t iterate_vectors_RI_Kmer(int max_size){
  auto container = std::vector<RI_Kmer>(max_size);
  matrix_t distances{24, 24};
  
  auto container_it = container.begin();
  auto container_end = container.end();
  int i = 0;
  while(container_it != container_end){
    distances(i%24,i%24) = (*container_it).next_char_prob[2] + 5.0;
    ++container_it;
    i++;
  }

  return distances;
}

int main(int argc, char *argv[]){
  int max_size = std::stoi(argv[2]);
  //auto array = iterate_kmers_bench_dev(max_size);
  matrix_t array;
  if(argv[1] == 'i'){
    array = iterate_vectors_int(max_size);
  } else {
    array = iterate_vectors_RI_Kmer(max_size);
  }
  for(int i = 0; i < 10; i++){
    if(i < 10) std::cout << array(i,0);
  }
}