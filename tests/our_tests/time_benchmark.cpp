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
#include "distances/dvstar.hpp"

using RI_Kmer = container::RI_Kmer; 

void prettyPrint(size_t insert_time_fst, size_t insert_time_snd, size_t find_time_fst, size_t find_time_snd, 
      size_t iterate_time, size_t dvstar_time, int items_fst, int items_snd, std::string container){
  int string_length = container.length() + 24; 
  std::cout << std::string(string_length, '-') << std::endl;
  std::cout << "|           " << container << "           |" << std::endl; 
  std::cout << std::string(string_length, '-') << std::endl; 
  // std::cout << std::string(string_length, ' ') << std::endl;
  std::cout << "    First Insert" << std::endl;  
  std::cout << "Total time : " << insert_time_fst / 1000 << " [micro sec] " << std::endl; 
  std::cout << "Sec / item : " << insert_time_fst / items_fst << " [nano sec] " << std::endl; 

  std::cout << "    Second Insert" << std::endl;  
  std::cout << "Total time : " << insert_time_snd / 1000 << " [micro sec] " << std::endl; 
  std::cout << "Sec / item : " << insert_time_snd / items_snd << " [nano sec] " << std::endl;

  std::cout << "    First Find" << std::endl;  
  std::cout << "Total time : " << find_time_fst / 1000 << " [micro sec] " << std::endl; 
  std::cout << "Sec / item : " << find_time_fst / items_fst << " [nano sec] " << std::endl;

  std::cout << "    Second Find" << std::endl;  
  std::cout << "Total time : " << find_time_snd / 1000 << " [micro sec] " << std::endl; 
  std::cout << "Sec / item : " << find_time_snd / items_snd << " [nano sec] " << std::endl;

  std::cout << "    Iterate" << std::endl;  
  std::cout << "Total time : " << iterate_time / 1000 << " [micro sec] " << std::endl; 
  std::cout << "Sec / item : " << iterate_time / (items_fst + items_snd) << " [nano sec] " << std::endl;

  std::cout << "    Dvstar" << std::endl;  
  std::cout << "Total time : " << dvstar_time / 1000 << " [micro sec] " << std::endl; 
  std::cout << "Sec / item : " << dvstar_time / (items_fst + items_snd) << " [nano sec] " << std::endl;
  std::cout << std::endl; 
}

template <typename VC>
void iterate_kmers_f(VC left, VC right){
  double dot_product = 0; 
  double left_norm = 0; 
  double right_norm = 0; 
  int background_order = 0; 

  left.iterate_kmers(left, right, [&](auto &left_v, auto &right_v) { 
    if (left_v.length <= background_order) {
      return;
    }
    auto [left_comp, right_comp] = distance::get_components(left_v, right_v);
    for (int i = 0; i < 4; i++) { 
      dot_product += left_comp[i] * right_comp[i];
      left_norm += std::pow(left_comp[i], 2.0);
      right_norm += std::pow(right_comp[i], 2.0);
    }
  });

  std::cout << "Dot product = " << dot_product << ", left norm = " << left_norm << ", right norm = " << right_norm << std::endl; 
}

std::tuple<std::vector<RI_Kmer>, int> get_kmer_vector(std::filesystem::path path){
  std::vector<RI_Kmer> kmers{};

  std::ifstream ifs(path, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);
  vlmc::VLMCKmer input_kmer{};

  int items = 0;

  while (ifs.peek() != EOF){
    archive(input_kmer);
    RI_Kmer ri_kmer(input_kmer); 
    kmers.push_back(ri_kmer);
    items++;
  }
  ifs.close();
  return std::make_tuple(kmers, items);
}

template <typename VC> 
void run_timer(std::string container){
  std::filesystem::path path_fst{"../data/one_human_VLMCs/human_genome_1.bintree"};
  std::filesystem::path path_snd{"../data/one_human_VLMCs/human_genome_2.bintree"};

  auto fst = get_kmer_vector(path_fst);
  auto snd = get_kmer_vector(path_snd);

  int items_fst = std::get<1>(fst);
  int items_snd = std::get<1>(snd); 
  auto kmers_fst = std::get<0>(fst);
  auto kmers_snd = std::get<0>(snd);

  std::chrono::steady_clock::time_point begin_insert_fst = std::chrono::steady_clock::now();
  VC vlmc_fst{path_fst}; 
  std::chrono::steady_clock::time_point end_insert_fst = std::chrono::steady_clock::now(); 
  auto insert_time_fst = std::chrono::duration_cast<std::chrono::nanoseconds>(end_insert_fst - begin_insert_fst).count();

  std::chrono::steady_clock::time_point begin_insert_snd = std::chrono::steady_clock::now();
  VC vlmc_snd{path_snd}; 
  std::chrono::steady_clock::time_point end_insert_snd = std::chrono::steady_clock::now(); 
  auto insert_time_snd = std::chrono::duration_cast<std::chrono::nanoseconds>(end_insert_snd - begin_insert_snd).count();

  auto begin_find_fst = std::chrono::steady_clock::now();
  for (int i = 0; i < kmers_fst.size(); i++){
    kmers_fst[i] = vlmc_fst.find(kmers_fst[i].integer_rep);
  }
  auto end_find_fst = std::chrono::steady_clock::now();
  auto find_time_fst = std::chrono::duration_cast<std::chrono::nanoseconds>(end_find_fst - begin_find_fst).count();

  auto begin_find_snd = std::chrono::steady_clock::now();
  for (int i = 0; i < kmers_snd.size(); i++){
    kmers_snd[i] = vlmc_snd.find(kmers_snd[i].integer_rep);
  }
  auto end_find_snd = std::chrono::steady_clock::now();
  auto find_time_snd = std::chrono::duration_cast<std::chrono::nanoseconds>(end_find_snd - begin_find_snd).count();

  auto begin_iterate = std::chrono::steady_clock::now();
  iterate_kmers_f(vlmc_fst, vlmc_snd);
  auto end_iterate = std::chrono::steady_clock::now();
  auto iterate_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_iterate - begin_iterate).count();

  auto begin_dvstar = std::chrono::steady_clock::now();
  distance::dvstar(vlmc_fst, vlmc_snd, 0.0);
  auto end_dvstar = std::chrono::steady_clock::now();
  auto dvstar_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_dvstar - begin_dvstar).count();

  auto items_fst_sum = 0; 
  auto items_snd_sum = 0; 
  for (int i = 0; i < kmers_fst.size(); i++){
    items_fst_sum += kmers_fst[i].integer_rep; 
  }
  for (int i = 0; i < kmers_snd.size(); i++){
    items_snd_sum += kmers_snd[i].integer_rep; 
  }

  std::cout << "Items in first = " << items_fst << " sums up to " << items_fst_sum << ", items in second = " << items_snd << " sums up to " << items_snd_sum << std::endl;
  
  prettyPrint(insert_time_fst, insert_time_snd, find_time_fst, find_time_snd, iterate_time, dvstar_time, items_fst, items_snd, container);
}

struct benchmark_kmer {
  size_t length = 0; 
  int integer_rep;
  // int background_rep; <- Should be implemented 
  std::array<double,4> next_char_prob{};
  std::array<uint64, 4> bit_representation; 
  bool is_null = true;
};

int get_index(const vlmc::VLMCKmer &kmer) {
  int integer_value = 0;
  int offset = 1;
  for (int i = kmer.length - 1; i >= 0; i--) {
    uchar row = i >> 5;
    uchar pos_in_row = i & 31;
    uchar n_shift_pos_to_end = (62 - pos_in_row * 2);
    auto kmer_2_bits = ((kmer.kmer_data[row] >> n_shift_pos_to_end) & 3) + 1;
    integer_value += (kmer_2_bits * offset);
    offset = offset << 2;
  }
  return integer_value;
}

std::tuple<int, int> read_in_kmer(vlmc::VLMCKmer old_kmer){
  benchmark_kmer bmer{};

  auto begin_basic = std::chrono::steady_clock::now();
  bmer.length = old_kmer.length;
  double child_count = std::accumulate(old_kmer.next_symbol_counts.begin(), old_kmer.next_symbol_counts.end(), pseudo_count_amount * 4);
  bmer.next_char_prob = {(double(old_kmer.next_symbol_counts[0]) + pseudo_count_amount) / child_count,
         (double(old_kmer.next_symbol_counts[1]) + pseudo_count_amount) / child_count,
         (double(old_kmer.next_symbol_counts[2]) + pseudo_count_amount) / child_count,
         (double(old_kmer.next_symbol_counts[3]) + pseudo_count_amount) / child_count};
  bmer.is_null = false;
  auto end_basic = std::chrono::steady_clock::now();
  auto basic_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_basic - begin_basic).count();

  auto begin_integer_rep = std::chrono::steady_clock::now();
  bmer.integer_rep = get_index(old_kmer);
  auto end_integer_rep = std::chrono::steady_clock::now();
  auto integer_rep_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_integer_rep - begin_integer_rep).count();

  return std::make_tuple(basic_time, integer_rep_time);
}

void benchmark_read_in_kmer(){
  int total_basic = 0; 
  int total_integer_rep = 0; 
  int count = 0; 

  std::filesystem::path path{"../data/one_human_VLMCs/human_genome_1.bintree"};
  std::ifstream ifs(path, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);
  Kmer kmer{};
  while (ifs.peek() != EOF){
    archive(kmer);
    auto res = read_in_kmer(kmer); 
    total_basic += std::get<0>(res);
    total_integer_rep += std::get<1>(res);
    count++; 
  }
  ifs.close();

  std::cout << std::endl; 
  std::cout << "Total time : " << (total_basic + total_integer_rep) / 1000 << " [micro sec]" << std::endl; 
  std::cout << "Basic time : " << total_basic / 1000 << " [micro sec], Avg : " << total_basic / count << " [nano sec]" << std::endl;  
  std::cout << "Integer_rep time : " << total_integer_rep / 1000 << " [micro sec], Avg : " << total_integer_rep / count << " [nano sec]" << std::endl;  
}

void benchmark_container_inv_sqrt(){
  std::filesystem::path path{"../data/one_human_VLMCs/human_genome_1.bintree"};

  auto begin_old = std::chrono::steady_clock::now();
  //old
  container::VLMC_sorted_vector(path, 1, false);

  auto end_old = std::chrono::steady_clock::now();
  auto time_old = std::chrono::duration_cast<std::chrono::microseconds>(end_old - begin_old).count();

  auto begin_new = std::chrono::steady_clock::now();
  //new
  container::VLMC_sorted_vector(path, 1, true);

  auto end_new = std::chrono::steady_clock::now();
  auto time_new = std::chrono::duration_cast<std::chrono::microseconds>(end_new - begin_new).count();
  
  std::cout << std::endl;
  std::string descriptive_string = "Time to construct sorted vector (eigen 1/sqrt)";
  int string_length = descriptive_string.length() + 24; 
  std::cout << std::string(string_length, '-') << std::endl;
  std::cout << "|           " << descriptive_string << "           |" << std::endl; 
  std::cout << std::string(string_length, '-') << std::endl; 
  std::cout << "Total time for old : " << time_old << " [micro sec]" << std::endl;  
  std::cout << "Total time for eigen : " << time_new << " [micro sec]" << std::endl; 
}

void benchmark_kmer_comparison(){
  int nr_kmers = 0;

  std::filesystem::path path{"../data/one_human_VLMCs/human_genome_1.bintree"};
  std::ifstream ifs(path, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);

  Kmer old_kmer{};

  auto begin_old_comp = std::chrono::steady_clock::now();

  std::vector<Kmer> old_kmers{};
  while (ifs.peek() != EOF){
    archive(old_kmer);
    old_kmers.push_back(old_kmer);
    nr_kmers++;
  }
  ifs.close();

  for(auto first_kmer : old_kmers){
    for(auto second_kmer : old_kmers){
      // Compare two old kmers;
      auto test_old = first_kmer < second_kmer;

    }
  }
  auto end_old_comp = std::chrono::steady_clock::now();
  auto old_kmer_comp_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_old_comp - begin_old_comp).count();

  auto begin_our_comp = std::chrono::steady_clock::now();

  container::VLMC_vector our_kmers{path};
  for(auto first_kmer : our_kmers){
    for(auto second_kmer : our_kmers){
      auto test_our = first_kmer < second_kmer;
    }
  }

  auto end_our_comp = std::chrono::steady_clock::now();
  auto our_kmer_comp_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_our_comp - begin_our_comp).count();
  
  std::cout << std::endl;
  std::string descriptive_string = "Time to compare < operator";
  int string_length = descriptive_string.length() + 24; 
  std::cout << std::string(string_length, '-') << std::endl;
  std::cout << "|           " << descriptive_string << "           |" << std::endl; 
  std::cout << std::string(string_length, '-') << std::endl; 
  std::cout << "Total time for old : " << old_kmer_comp_time << " [nano sec], Avg : " << old_kmer_comp_time / nr_kmers << " [nano sec]" << std::endl;  
  std::cout << "Total time for our : " << our_kmer_comp_time << " [nano sec], Avg : " << our_kmer_comp_time / nr_kmers << " [nano sec]" << std::endl;  
}

int main(int argc, char *argv[]){
  // int num_items = 1500;

  // run_timer<container::VLMC_vector>("Vector");
  run_timer<container::VLMC_Indexing>("Indexing");
  run_timer<container::VLMC_sorted_vector>("Sorted Vector");
  run_timer<container::VLMC_B_tree>("B-tree");
  run_timer<container::VLMC_hashmap>("Hashmap");
  run_timer<container::VLMC_Combo>("Combo");
  run_timer<container::VLMC_Veb>("Veb-tree");
  benchmark_read_in_kmer();
  benchmark_kmer_comparison();
  benchmark_container_inv_sqrt();
}