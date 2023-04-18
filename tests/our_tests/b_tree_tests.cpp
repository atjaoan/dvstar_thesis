#pragma once 

#include <gtest/gtest.h>

#include <functional>
#include <filesystem>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <stdlib.h> 

#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_container.hpp"
#include "read_in_kmer.hpp"
#include "global_aliases.hpp"
#include "b_tree.hpp"
#include "b_tree_layout.hpp"

using RI_Kmer = container::RI_Kmer; 

class BTreeTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path path_to_bintrees{"../data/test_VLMCs"};
};

b_tree::B_Tree createTree(std::filesystem::path path){
  int background_order = 0;
  eigenx_t cached_context((int)std::pow(4, background_order), 4);

  std::vector<RI_Kmer> container{};

  auto fun = [&](const RI_Kmer &kmer) { container.push_back(kmer); }; 

  int offset_to_remove = container::load_VLMCs_from_file(path, cached_context, fun, background_order);
      
  std::sort(std::execution::seq, container.begin(), container.end());
  for (size_t i = 0; i < container.size(); i++){
    RI_Kmer kmer = container[i];
    int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
    int offset = background_idx - offset_to_remove;
    container[i].next_char_prob *= cached_context.row(offset).rsqrt();
  }
  std::cout << "Size " << container.size() << std::endl; 
  return b_tree::B_Tree(container);
}

b_tree::B_Tree createDummyTree(int x){
  std::vector<RI_Kmer> tmp{};

  for (int i = 0; i < x; i++){
    tmp.push_back(RI_Kmer(i));
  }

  return b_tree::B_Tree(tmp);
}

/*
std::tuple<std::vector<RI_Kmer>, std::vector<RI_Kmer>> generateRandomVector(int max_size, int match_prc) {
  std::vector<RI_Kmer> left;
  std::vector<RI_Kmer> right; 
  int count = 0; 
  int matches = 0; 
  int number = 0; 
  bool last_left = false;
  while (count < max_size) {
    RI_Kmer new_kmer = RI_Kmer(number, {{0.25, 0.25, 0.25, 0.25}});
    int match_chance = rand() % 100;
    if (match_chance < match_prc) {
      left.push_back(new_kmer);
      right.push_back(new_kmer);
      count++; 
      matches++; 
    } else {
      if (last_left) {
        right.push_back(new_kmer);
        last_left = false;
        count++;  
      } else {
        left.push_back(new_kmer);
        last_left = true; 
      }
    }
    number++; 
  }
  std::cout << "The count is = " << count << " and matches = " << matches << " prc = " << float(matches) / float(count) << std::endl; 
  return std::make_tuple(left, right); 
}
*/

const int B = 3; 

float measure_distance(b_tree::B_Tree t1, b_tree::B_Tree t2) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;

  int left_nblocks = t1.get_nblocks();

  int left_k = 0;
  int left_i = 0; 

  int count = 0;
  int count_second_hit = 0;  

  while (left_k < left_nblocks){
    RI_Kmer &left_kmer = t1.get(left_k, left_i);
    int search_k = std::get<0>(t1.go_back(left_k));
    auto right_res = t2.search(left_kmer.integer_rep, search_k); 

    if (std::get<0>(right_res) && left_kmer == std::get<1>(right_res)){
      dot_product += (left_kmer.next_char_prob * std::get<1>(right_res).next_char_prob).sum();
      left_norm += left_kmer.next_char_prob.square().sum();
      right_norm += std::get<1>(right_res).next_char_prob.square().sum();
      count++; 
    } else {
      auto right_res = t2.search(left_kmer.integer_rep, 0);
      RI_Kmer &right_kmer = std::get<1>(right_res);

      if (std::get<0>(right_res) && left_kmer == std::get<1>(right_res)){
        dot_product += (left_kmer.next_char_prob * std::get<1>(right_res).next_char_prob).sum();
        left_norm += left_kmer.next_char_prob.square().sum();
        right_norm += std::get<1>(right_res).next_char_prob.square().sum();
        count_second_hit++; 
        count++;
      } 
    }
    if (left_i < B - 1) {
      left_i++;
    } else {
      left_k++;
      left_i = 0; 
    }
  }

  std::cout << "Total hits = " << count << " second hit = " << count_second_hit << std::endl; 

  return container::normalise_dvstar(dot_product, left_norm, right_norm);
}

int go(int k, int i) { return k * (B + 1) + i + 1; }

std::tuple<int,int> go_back(int k) { 
  int old_i = (k - 1) % (B + 1); 
  int old_k = (k - old_i - 1) / (B + 1);
  return std::make_tuple(old_k, old_i);  
}

void increment(int& k, int& i, int& level, const int nblocks) {
  if (level > 0){
    k = go(k, i + 1);
    i = 0;
    level--; 
    while(level > 0){
      k = go(k, i); 
      level--;
    }
    if (k > nblocks - 1) {
      auto res = go_back(k);
      k = std::get<0>(res);
      i = std::get<1>(res);
    }
  } else if (i < B - 1){
    i++; 
  } else {
    auto res = go_back(k);
    k = std::get<0>(res);
    i = std::get<1>(res);  
    level++;
    if (i == B){
      while (i == B){
        auto tmp = go_back(k);
        k = std::get<0>(tmp);
        i = std::get<1>(tmp);
        level++;
      }
    }
  }
}

float crazy_search_iteration(b_tree::B_Tree left_kmers, b_tree::B_Tree right_kmers) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;

  auto left_res = left_kmers.search(1);
  auto right_res = right_kmers.search(1); 

  int left_nblocks = left_kmers.get_nblocks();
  int right_nblocks = right_kmers.get_nblocks();

  int left_k = std::get<0>(std::get<1>(left_res));
  int right_k = std::get<0>(std::get<1>(right_res)); 

  int left_i = std::get<1>(std::get<1>(left_res));
  int right_i = std::get<1>(std::get<1>(right_res)); 

  int left_level = 0;
  int right_level = 0; 

  int count = 0; 
  int count_hit = 0; 

  while(left_i >= 0 && left_k < left_nblocks && right_i >= 0 && right_k < right_nblocks){
    RI_Kmer &left_kmer = left_kmers.get(left_k, left_i);
    RI_Kmer &right_kmer = right_kmers.get(right_k, right_i);

    if(left_kmer == right_kmer){
      dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
      left_norm += left_kmer.next_char_prob.square().sum();
      right_norm += right_kmer.next_char_prob.square().sum();
      increment(left_k, left_i, left_level, left_nblocks);
      increment(right_k, right_i, right_level, right_nblocks);
      count++; 
    } else if (left_kmer > right_kmer) {
      auto right_res = right_kmers.search(left_kmer.integer_rep);
      if (std::get<0>(right_res) > -1){
        right_level = std::get<0>(right_res);
        right_k = std::get<0>(std::get<1>(right_res));
        right_i = std::get<1>(std::get<1>(right_res));
        RI_Kmer &right_kmer = right_kmers.get(right_k, right_i);
        if (left_kmer == right_kmer){
          dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
          left_norm += left_kmer.next_char_prob.square().sum();
          right_norm += right_kmer.next_char_prob.square().sum();
          increment(left_k, left_i, left_level, left_nblocks);
          increment(right_k, right_i, right_level, right_nblocks);
          count++; 
          count_hit++;
        }
      } else {
        increment(left_k, left_i, left_level, left_nblocks);
        increment(right_k, right_i, right_level, right_nblocks);
      }
    } else {
      auto left_res = left_kmers.search(right_kmer.integer_rep);
      if (std::get<0>(left_res) > -1){
        left_level = std::get<0>(left_res);
        left_k = std::get<0>(std::get<1>(left_res));
        left_i = std::get<1>(std::get<1>(left_res));
        RI_Kmer &left_kmer = left_kmers.get(left_k, left_i);
        if (left_kmer == right_kmer){
          dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
          left_norm += left_kmer.next_char_prob.square().sum();
          right_norm += right_kmer.next_char_prob.square().sum();
          increment(left_k, left_i, left_level, left_nblocks);
          increment(right_k, right_i, right_level, right_nblocks);
          count++; 
          count_hit++;
        }
      } else {
        increment(left_k, left_i, left_level, left_nblocks);
        increment(right_k, right_i, right_level, right_nblocks);
      }
    }
  }

  std::cout << "The count is " << count << std::endl; 
  std::cout << "The count for search hits is = " << count_hit << std::endl; 

  return container::normalise_dvstar(dot_product, left_norm, right_norm);
}

float crazy_measure_distance(b_tree::B_Tree left_kmers, b_tree::B_Tree right_kmers) {
  out_t dot_product = 0.0;
  out_t left_norm = 0.0;
  out_t right_norm = 0.0;

  auto left_res = left_kmers.search(1);
  auto right_res = right_kmers.search(1); 

  int left_nblocks = left_kmers.get_nblocks();
  int right_nblocks = right_kmers.get_nblocks();

  int left_k = std::get<0>(std::get<1>(left_res));
  int right_k = std::get<0>(std::get<1>(right_res)); 

  int left_i = std::get<1>(std::get<1>(left_res));
  int right_i = std::get<1>(std::get<1>(right_res)); 

  int left_level = 0;
  int right_level = 0; 

  int count = 0; 
  int count_hit = 0; 

  int max_level_left = 0; 
  int max_level_right = 0; 

  while(left_i >= 0 && left_k < left_nblocks && right_i >= 0 && right_k < right_nblocks){
    if (max_level_left < left_level){
      max_level_left = left_level;
    }
    if (max_level_right < right_level){
      max_level_right = right_level; 
    }

    RI_Kmer &left_kmer = left_kmers.get(left_k, left_i);
    RI_Kmer &right_kmer = right_kmers.get(right_k, right_i);

    if(left_kmer == right_kmer){
      dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
      left_norm += left_kmer.next_char_prob.square().sum();
      right_norm += right_kmer.next_char_prob.square().sum();
      increment(left_k, left_i, left_level, left_nblocks);
      increment(right_k, right_i, right_level, right_nblocks);
      count++; 
    } else if (left_kmer > right_kmer) {
      increment(right_k, right_i, right_level, right_nblocks);
      if (right_i >= 0 && right_k < right_nblocks){
        RI_Kmer &right_kmer = right_kmers.get(right_k, right_i);
        if (left_kmer == right_kmer){
          dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
          left_norm += left_kmer.next_char_prob.square().sum();
          right_norm += right_kmer.next_char_prob.square().sum();
          increment(left_k, left_i, left_level, left_nblocks);
          increment(right_k, right_i, right_level, right_nblocks);
          count++; 
        } else {
          auto right_res = right_kmers.search(left_kmer.integer_rep);
          if (std::get<0>(right_res) > -1){
            right_level = std::get<0>(right_res);
            right_k = std::get<0>(std::get<1>(right_res));
            right_i = std::get<1>(std::get<1>(right_res));
            RI_Kmer &right_kmer = right_kmers.get(right_k, right_i);
            if (left_kmer == right_kmer){
              dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
              left_norm += left_kmer.next_char_prob.square().sum();
              right_norm += right_kmer.next_char_prob.square().sum();
              increment(left_k, left_i, left_level, left_nblocks);
              increment(right_k, right_i, right_level, right_nblocks);
              count++; 
              count_hit++;
            }
          } else {
            increment(left_k, left_i, left_level, left_nblocks);
            increment(right_k, right_i, right_level, right_nblocks);
          }
        }
      }
    } else {
      increment(left_k, left_i, left_level, left_nblocks);
      if (left_i >= 0 && left_k < left_nblocks){
        RI_Kmer &left_kmer = left_kmers.get(left_k, left_i);
        if (left_kmer == right_kmer){
          dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
          left_norm += left_kmer.next_char_prob.square().sum();
          right_norm += right_kmer.next_char_prob.square().sum();
          increment(left_k, left_i, left_level, left_nblocks);
          increment(right_k, right_i, right_level, right_nblocks);
          count++; 
        } else {
          auto left_res = left_kmers.search(right_kmer.integer_rep);
          if (std::get<0>(left_res) > -1){
            left_level = std::get<0>(left_res);
            left_k = std::get<0>(std::get<1>(left_res));
            left_i = std::get<1>(std::get<1>(left_res));
            RI_Kmer &left_kmer = left_kmers.get(left_k, left_i);
            if (left_kmer == right_kmer){
              dot_product += (left_kmer.next_char_prob * right_kmer.next_char_prob).sum();
              left_norm += left_kmer.next_char_prob.square().sum();
              right_norm += right_kmer.next_char_prob.square().sum();
              increment(left_k, left_i, left_level, left_nblocks);
              increment(right_k, right_i, right_level, right_nblocks);
              count++; 
              count_hit++;
            }
          } else {
            increment(left_k, left_i, left_level, left_nblocks);
            increment(right_k, right_i, right_level, right_nblocks);
          }
        }
      }
    }
  }

  std::cout << "Max level left = " << max_level_left << " nblocks = " << left_nblocks << std::endl; 
  std::cout << "Max level right = " << max_level_right << " nblocks = " << right_nblocks << std::endl; 

  std::cout << "The count is " << count << std::endl; 
  std::cout << "The count for search hits is = " << count_hit << std::endl; 

  return container::normalise_dvstar(dot_product, left_norm, right_norm);
}

TEST_F(BTreeTests, increment_test) {
  int size = 5456;
  b_tree::B_Tree tree = createDummyTree(size);

  const int nblocks = tree.get_nblocks();

  int level = 0;
  int k = std::get<0>(std::get<1>(tree.search(0)));
  int i = std::get<1>(std::get<1>(tree.search(0))); 

  for (int x = 0; x < size; x++){
    int res = tree.get(k, i).integer_rep;
    //std::cout << "X is " << x << " vs " << tree.get(k, i).integer_rep << std::endl; 
    EXPECT_EQ(x, res);
    increment(k, i, level, nblocks);
  }
}

TEST_F(BTreeTests, new_btree) {
  std::filesystem::path path_bintree_1{"../data/benchmarking/human/large/human_genome_1_2_3_10.bintree"};
  std::filesystem::path path_bintree_2{"../data/benchmarking/human/large/human_genome_10_2_3_10.bintree"};

  // std::chrono::steady_clock::time_point start_create_tree = std::chrono::steady_clock::now();
  b_tree::B_Tree tree1 = createTree(path_bintree_1);
  b_tree::B_Tree tree2 = createTree(path_bintree_2);
  // std::chrono::steady_clock::time_point stop_create_tree = std::chrono::steady_clock::now();
// 
  // std::chrono::steady_clock::time_point start_create_container = std::chrono::steady_clock::now();
  container::VLMC_sorted_vector vec1 = container::VLMC_sorted_vector(path_bintree_1, 0);
  container::VLMC_sorted_vector vec2 = container::VLMC_sorted_vector(path_bintree_2, 0);
  // std::chrono::steady_clock::time_point stop_create_container = std::chrono::steady_clock::now();

  // auto create_tree = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_create_tree - start_create_tree).count();
  // auto create_container = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_create_container - start_create_container).count();
// 
  // std::cout << "Time in nanoseconds create tree      = " << create_tree << std::endl; 
  // std::cout << "Time in nanoseconds create container = " << create_container << std::endl; 

  // auto generated_vecs = generateRandomVector(1000000, 2);

  std::cout << "Generate vectors" << std::endl;  
  // b_tree::B_Tree tree1 = b_tree::B_Tree(std::get<0>(generated_vecs));
  // b_tree::B_Tree tree2 = b_tree::B_Tree(std::get<1>(generated_vecs));
  // container::VLMC_sorted_vector vec1 = container::VLMC_sorted_vector(std::get<0>(generated_vecs));
  // container::VLMC_sorted_vector vec2 = container::VLMC_sorted_vector(std::get<1>(generated_vecs));

  std::cout << "Finished creation" << std::endl; 

  // b_tree::B_Tree tree1 = createDummyTree();

  std::chrono::steady_clock::time_point start_measure_distance = std::chrono::steady_clock::now();
  auto dist = measure_distance(tree1, tree2);
  std::chrono::steady_clock::time_point stop_measure_distance = std::chrono::steady_clock::now();

  std::chrono::steady_clock::time_point start_crazy_measure = std::chrono::steady_clock::now();
  auto c_dist = crazy_measure_distance(tree1, tree2);
  std::chrono::steady_clock::time_point stop_crazy_measure = std::chrono::steady_clock::now();

  std::chrono::steady_clock::time_point start_crazy_search = std::chrono::steady_clock::now();
  auto s_dist = crazy_search_iteration(tree1, tree2);
  std::chrono::steady_clock::time_point stop_crazy_search = std::chrono::steady_clock::now();

  std::chrono::steady_clock::time_point start_iterate_kmers = std::chrono::steady_clock::now();
  auto v_dist = container::iterate_kmers(vec1, vec2);
  std::chrono::steady_clock::time_point stop_iterate_kmers = std::chrono::steady_clock::now();

  auto measure_distance = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_measure_distance - start_measure_distance).count();
  auto crazy_measure_distance = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_crazy_measure - start_crazy_measure).count();
  auto iterate_kmers = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_iterate_kmers - start_iterate_kmers).count();
  auto crazy_search = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_crazy_search - start_crazy_search).count();

  std::cout << "Distance " << dist << " vs " << c_dist << " vs " << v_dist << " vs " << s_dist << std::endl; 
  std::cout << "Time in milliseconds measure distance = " << measure_distance / 1000 << std::endl; 
  std::cout << "Time in milliseconds crazy measure    = " << crazy_measure_distance / 1000 << std::endl; 
  std::cout << "Time in milliseconds iterate kmers    = " << iterate_kmers / 1000 << std::endl; 
  std::cout << "Time in milliseconds crazy search     = " << crazy_search / 1000 << std::endl; 
}
