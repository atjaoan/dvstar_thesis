#pragma once 

#include <functional>
#include <filesystem>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_container.hpp"
#include "cluster_container.hpp"
#include "../read_helper.hpp"
#include "read_in_kmer.hpp"

void prettyPrint(size_t insert_time, size_t find_time, size_t iterate_time, int items, std::string container, int control_count){
  int string_length = container.length() + 24; 
  std::cout << std::string(string_length, '-') << std::endl;
  std::cout << "|           " << container << "           |" << std::endl; 
  std::cout << std::string(string_length, '-') << std::endl; 
  std::cout << std::string(string_length, ' ') << std::endl;
  std::cout << "    Insert" << std::endl;  
  std::cout << "Total time : " << insert_time << " [nano sec] " << std::endl; 
  std::cout << "Sec / item : " << insert_time  / items << " [nano sec] " << std::endl; 
  std::cout << std::string(string_length, ' ') << std::endl;

  std::cout << "    Find" << std::endl;  
  std::cout << "Total time : " << find_time << " [nano sec] " << std::endl; 
  std::cout << "Sec / item : " << find_time  / items << " [nano sec] " << std::endl;
  std::cout << std::string(string_length, ' ') << std::endl;

  std::cout << "    Iterate" << std::endl;  
  std::cout << "Total time : " << iterate_time << " [nano sec] " << std::endl; 
  std::cout << "Sec / item : " << iterate_time  / items << " [nano sec] " << std::endl;
  std::cout << "Control count (" << control_count << ")" << std::endl; 
}

template <typename VC> 
void run_timer(int items, VC vlmc, std::string container){
  std::chrono::steady_clock::time_point begin_insert = std::chrono::steady_clock::now();
  for (size_t i = 0; i < items; i++){
    auto kmer = container::RI_Kmer{i};
    vlmc.push(kmer);
  }
  std::chrono::steady_clock::time_point end_insert = std::chrono::steady_clock::now(); 
  auto insert_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_insert - begin_insert).count();

  auto begin_find = std::chrono::steady_clock::now();
  for (size_t i = 0; i < items; i++){
    vlmc.find(i);
  }
  auto end_find = std::chrono::steady_clock::now();
  auto find_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_find - begin_find).count();

  int control_count = 0; 

  auto begin_iterate = std::chrono::steady_clock::now();
  for (size_t i = 0; i < items; i++){
    vlmc.iterate_kmers(vlmc, vlmc, [&](auto &, auto &) { control_count++; });
  }
  auto end_iterate = std::chrono::steady_clock::now();
  auto iterate_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_iterate - begin_iterate).count();
  
  prettyPrint(insert_time, find_time, iterate_time, items, container, control_count);
}

int main(int argc, char *argv[]){
  int num_items = 1500;

  run_timer(num_items, container::VLMC_vector{}, "Vector");

  run_timer(num_items, container::Index_by_value{}, "Index by Value");

  run_timer(num_items, container::VLMC_sorted_vector{}, "Sorted Vector");

  run_timer(num_items, container::VLMC_B_tree{}, "B-tree");

  run_timer(num_items, container::VLMC_hashmap{}, "Hashmap");
}