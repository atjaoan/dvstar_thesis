#pragma once 

#include <functional>
#include <filesystem>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>

#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_container.hpp"
#include "cluster_container.hpp"
#include "../read_helper.hpp"
#include "read_in_kmer.hpp"
#include "distances/dvstar.hpp"
#include "calc_dists.hpp"
#include "get_cluster.hpp"

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;
using RI_Kmer = container::RI_Kmer; 
using matrix_t = Eigen::MatrixXf;

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
  std::filesystem::path path_fst{"../data/test_VLMCs/sequences_1.bintree"};
  std::filesystem::path path_snd{"../data/test_VLMCs/sequences_2.bintree"};

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

void benchmark_kmer_major_load_calc(int nr_cores){
  std::cout << "Running test with " << nr_cores << " cores" << std::endl; 
  std::filesystem::path path{"../data/medium_test"};
  auto start_loading = std::chrono::steady_clock::now();
  auto cluster = cluster::get_kmer_cluster(path, 0);
  auto end_loading = std::chrono::steady_clock::now();

  matrix_t distance_matrix = calculate::calculate_distance_major(cluster, cluster, nr_cores);
  auto end_calc = std::chrono::steady_clock::now();

  auto time_loading = std::chrono::duration_cast<std::chrono::microseconds>(end_loading - start_loading).count();
  auto time_calc = std::chrono::duration_cast<std::chrono::microseconds>(end_calc - end_loading).count();
  auto time_total = std::chrono::duration_cast<std::chrono::microseconds>(end_calc - start_loading).count();
  
  std::cout << std::endl;
  std::string descriptive_string = "Time calc kmer-major: ";
  int string_length = descriptive_string.length() + 24; 
  std::cout << std::string(string_length, '-') << std::endl;
  std::cout << "|           " << descriptive_string << "           |" << std::endl; 
  std::cout << std::string(string_length, '-') << std::endl; 
  std::cout << "Time loading : " << time_loading << " [micro sec]" << " , fraction: " << (time_loading / (double)time_total) << std::endl;  
  std::cout << "Time calc : " << time_calc << " [micro sec]"  << " , fraction: " << (time_calc / (double)time_total) << std::endl; 
  std::cout << "Total time : " << time_total << " [micro sec]" << std::endl; 
}

/*
void calculate_kmer_buckets(size_t start_bucket, size_t stop_bucket, 
    matrix_t &distances, matrix_t &dot_prod, matrix_t &left_norm, matrix_t &right_norm,
    container::Kmer_Cluster &cluster_left, container::Kmer_Cluster &cluster_right) {
  int get = 0;
  int dvstar = 0;
  int total = 0;
  for (size_t i = start_bucket; i < stop_bucket; i++) {
    if (cluster_left.is_bucket_empty(i) || cluster_right.is_bucket_empty(i)){
      continue; 
    }
    auto start_get = std::chrono::steady_clock::now();
    auto vec_left_begin = cluster_left.get_bucket_begin(i);
    auto vec_left_end = cluster_left.get_bucket_end(i);
    auto vec_right_begin = cluster_right.get_bucket_begin(i);
    auto vec_right_end = cluster_right.get_bucket_end(i);
    auto end_get = std::chrono::steady_clock::now();

    distance::dvstar_kmer_major(vec_left_begin, vec_left_end, vec_right_begin, vec_right_end, dot_prod, left_norm, right_norm);
    auto end_dist = std::chrono::steady_clock::now();

    auto time_get = std::chrono::duration_cast<std::chrono::microseconds>(end_get - start_get).count();
    auto time_dvstar = std::chrono::duration_cast<std::chrono::microseconds>(end_dist - end_get).count();
    auto time_both = std::chrono::duration_cast<std::chrono::microseconds>(end_dist - start_get).count();

    get += time_get;
    dvstar += time_dvstar;
    total += time_both;
  }
  std::cout << "Time get : " << get << " [micro sec]" << " , fraction: " << (get / (double)total) << std::endl;  
  std::cout << "Time dvstar : " << dvstar << " [micro sec]"  << " , fraction: " << (dvstar / (double)total) << std::endl; 
  std::cout << "Total time in buckets : " << total << " [micro sec]" << std::endl; 
  std::cout << "I work on buckets : " << start_bucket << " to " << stop_bucket << std::endl;
}
*/


void benchmark_calculate_distance_major(){
  std::filesystem::path path{"../data/test_VLMCs"};
  int requested_cores = 4;
  auto cluster_left = cluster::get_kmer_cluster(path, 0);
  auto cluster_right = cluster::get_kmer_cluster(path, 0);

  const size_t processor_count = std::thread::hardware_concurrency();
  size_t used_cores = 1;
  if(requested_cores > cluster_left.experimental_bucket_count()){
    used_cores = cluster_left.experimental_bucket_count();
  } else if(requested_cores <= processor_count){
      used_cores = requested_cores;
  } else {
    used_cores = processor_count;
  }
  BS::thread_pool pool(used_cores);

  auto start_create_matricies = std::chrono::steady_clock::now();
  matrix_t distances = matrix_t::Zero(cluster_left.size(), cluster_right.size());
  matrix_t dot_prod = matrix_t::Zero(cluster_left.size(), cluster_right.size());
  matrix_t left_norm = matrix_t::Zero(cluster_left.size(), cluster_right.size());
  matrix_t right_norm = matrix_t::Zero(cluster_left.size(), cluster_right.size());
  auto end_create_matricies = std::chrono::steady_clock::now();

  auto start_dist = std::chrono::steady_clock::now();
  
  std::vector<matrix_t> dot_prods{};
  std::vector<matrix_t> lnorms{};
  std::vector<matrix_t> rnorms{};

  auto fun = [&](size_t start_bucket, size_t stop_bucket) {
    matrix_t dot_local = matrix_t::Zero(cluster_left.size(), cluster_right.size());
    matrix_t lnorm_local = matrix_t::Zero(cluster_left.size(), cluster_right.size());
    matrix_t rnorm_local = matrix_t::Zero(cluster_left.size(), cluster_right.size());
    calculate::calculate_kmer_buckets(start_bucket, stop_bucket, dot_local, lnorm_local, rnorm_local, cluster_left, cluster_right);
    dot_prods.push_back(dot_local);
    lnorms.push_back(lnorm_local);
    rnorms.push_back(rnorm_local);
  };
  
  parallel::pool_parallelize(cluster_left.experimental_bucket_count(), fun, requested_cores, pool);

  auto start_add = std::chrono::steady_clock::now();

  for(int i = 0; i < dot_prods.size(); ++i){
    dot_prod += dot_prods[i];
    left_norm += lnorms[i];
    right_norm += rnorms[i];
  }
  
  auto end_dist = std::chrono::steady_clock::now();

  auto start_norm = std::chrono::steady_clock::now();
  
  auto rec_fun = [&](size_t left, size_t right) {
    distances(left, right) = distance::normalise_dvstar(dot_prod(left, right), left_norm(left, right), right_norm(left, right));
  }; 

  auto norm_fun = [&](size_t start_vlmc, size_t stop_vlmc) {
    utils::matrix_recursion(start_vlmc, stop_vlmc, 0, cluster_right.size(), rec_fun); 
  };

  parallel::parallelize(cluster_left.size(), norm_fun, 1); 

  parallel::parallelize(cluster_left.size(), norm_fun, requested_cores);
  auto end_norm = std::chrono::steady_clock::now();

  auto time_create = std::chrono::duration_cast<std::chrono::microseconds>(end_create_matricies - start_create_matricies).count();
  auto time_dist = std::chrono::duration_cast<std::chrono::microseconds>(end_dist - start_dist).count();
  auto time_add = std::chrono::duration_cast<std::chrono::microseconds>(end_dist - start_add).count();
  auto time_norm = std::chrono::duration_cast<std::chrono::microseconds>(end_norm - start_norm).count();
  auto time_tot = std::chrono::duration_cast<std::chrono::microseconds>(end_norm - start_create_matricies).count();
  
  std::cout << std::endl;
  std::string descriptive_string = "Time calculate_distance_major: ";
  int string_length = descriptive_string.length() + 24; 
  std::cout << std::string(string_length, '-') << std::endl;
  std::cout << "|           " << descriptive_string << "           |" << std::endl; 
  std::cout << std::string(string_length, '-') << std::endl; 
  std::cout << "Time create matricies : " << time_create << " [micro sec]" << " , fraction: " << (time_create / (double)time_tot) << std::endl;  
  std::cout << "Time calculate dist : " << time_dist << " [micro sec]"  << " , fraction: " << (time_dist / (double)time_tot) << std::endl; 
  std::cout << "Time add matricies : " << time_add << " [micro sec]"  << " , fraction: " << (time_add / (double)time_dist) << std::endl; 
  std::cout << "Time normalize dist : " << time_norm << " [micro sec]"  << " , fraction: " << (time_norm / (double)time_tot) << std::endl; 
  std::cout << "Total time : " << time_tot << " [micro sec]" << std::endl; 
}

void print_kmers_to_file(){
  std::filesystem::path path_vlmcs{"../data/one_vlmcs"};
  std::filesystem::path output_path{"../tmp/" + path_vlmcs.filename().string() + "_kmer-distribution.txt"};
  utils::output_kmer_reps_to_file(path_vlmcs, output_path);
}

void get_kmer_vector_test(std::filesystem::path path, std::vector<RI_Kmer>& kmers){
  std::ifstream ifs(path, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);
  vlmc::VLMCKmer input_kmer{};

  while (ifs.peek() != EOF){
    archive(input_kmer);
    RI_Kmer ri_kmer(input_kmer); 
    kmers.push_back(ri_kmer);
  }
  ifs.close();
}

constexpr int reps = 20;
std::array<double, reps> iterate_kmers_bench_old(){
  std::filesystem::path path_fst{"../data/human_VLMCs/human_genome_1.bintree"};
  std::filesystem::path path_snd{"../data/human_VLMCs/human_genome_2.bintree"};

  std::vector<RI_Kmer> left{};
  std::vector<RI_Kmer> right{};
  get_kmer_vector_test(path_fst, left);
  get_kmer_vector_test(path_snd, right);

  std::sort(left.begin(), left.end());
  std::sort(right.begin(), right.end());
  
  std::array<double, reps> dot_products;
  auto start = std::chrono::steady_clock::now();
  for(int i = 0; i < reps; i++){
    double dot_product = 0.0;

    double left_norm = 0.0;
    double right_norm = 0.0;

    auto dvstar_fun = [&](auto &left_v, auto &right_v) {
      dot_product += (left_v.next_char_prob * right_v.next_char_prob).sum();
      left_norm += left_v.next_char_prob.square().sum();
      right_norm += right_v.next_char_prob.square().sum();
    };
    /*
    auto right_it = right.begin();
    auto right_end = right.end();
    for(auto &left_kmer : left){
      while((*right_it) < left_kmer){
        ++right_it;
        if(right_it == right_end) break;
      }
      if(left_kmer == (*right_it)){
        dvstar_fun(left_kmer, *right_it);
        ++right_it;
      }
    }
    */
    int left_i = 0;
    int right_i = 0;  
    const int left_size = left.size();
    const int right_size = right.size();
    while(left_i < left_size && right_i < right_size) {
      const RI_Kmer &left_kmer = left[left_i];
      const RI_Kmer &right_kmer = right[right_i];
      if (right_kmer == left_kmer) {
        dvstar_fun(left_kmer, right_kmer);
        left_i++;
        right_i++;
      } else if (left_kmer < right_kmer) {
        left_i++;
      } else {
        right_i++;
      }
    }
    
    dot_products[i] = dot_product;
  }
  auto end = std::chrono::steady_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  std::cout << "Time : " << time / (float) reps << " [micro sec]" << std::endl; 
  std::cout << "Dot products for sanity check: " << dot_products[0] << "\n";
  std::cout << "Size of RI_Kmer: " << sizeof(RI_Kmer) << " bytes" << "\n";
  std::cout << "Size of RI_Kmer in vector: " << sizeof(left[0]) << " bytes" << "\n";
  std::cout << "Size of left: " << sizeof(left) + sizeof(RI_Kmer) * left.capacity() << " bytes" << "\n";
  std::cout << "Size of right: " << sizeof(right) + sizeof(RI_Kmer) * right.capacity() << " bytes" << "\n";
  std::cout << "Address first element left: " << &left[0] << "\n";
  std::cout << "Address first element right: " << &right[0] << "\n";
  std::cout << "Address left: " << &left << "\n";
  std::cout << "Address right: " << &right << "\n";
  return dot_products;
}

// A function to implement bubble sort
void bubbleSort(std::vector<RI_Kmer> &container) {
  auto it_end = container.rend(); 
  for (int i = 0; i < container.size(); i++){
    for(auto it1 = container.rbegin(), it2 = ++container.rbegin(); it2 != it_end; it1++, it2++){
      if (*it1 > *it2){
        auto tmp = *it1;
        *it1 = *it2;
        *it2 = tmp;   
        // std::iter_swap(it1, it2);
      } 
    }
    it_end--; 
  }
}

std::tuple<int, int> run_sorting_timer(std::vector<std::vector<RI_Kmer>> containers, std::function<void(std::vector<RI_Kmer> &container)> sort_fun){
  int total_time = 0; 
  int avg_kmer_time = 0; 
  for (auto container : containers){
    std::vector<RI_Kmer> new_container(container);
    auto start_sort = std::chrono::steady_clock::now();
    sort_fun(new_container);
    auto end_sort = std::chrono::steady_clock::now();
    total_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end_sort - start_sort).count();
    avg_kmer_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end_sort - start_sort).count() / container.size();
  }
  return std::make_tuple(total_time / containers.size(), avg_kmer_time / containers.size()); 
}

void time_sorting_algorithms(std::filesystem::path directory){
  std::vector<std::tuple<std::tuple<int,int>,std::string>> sort_timings{};
  std::vector<std::vector<RI_Kmer>> containers{}; 

  auto sort_fun_bubblesort = [](std::vector<RI_Kmer> &container) {
    bubbleSort(container);
  };

  // auto sort_fun_timsort = [](std::vector<RI_Kmer> &container) {
  //   gfx::timsort(container.rbegin(), container.rend());
  // };

  auto sort_fun_seq = [](std::vector<RI_Kmer> &container) {
    std::sort(std::execution::seq, container.rbegin(), container.rend());
  };

  auto sort_fun_seq_not_rev = [](std::vector<RI_Kmer> &container) {
    std::sort(std::execution::seq, container.begin(), container.end());
  };

  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    std::vector<RI_Kmer> container{}; 
    std::ifstream ifs(dir_entry.path(), std::ios::binary);
    cereal::BinaryInputArchive archive(ifs);
    Kmer input_kmer{};
    while (ifs.peek() != EOF){
      archive(input_kmer);
      RI_Kmer ri_kmer{input_kmer};
      container.push_back(ri_kmer);
    }
    ifs.close();
    containers.push_back(container);  
  }
  sort_timings.push_back(std::make_tuple(run_sorting_timer(containers, sort_fun_seq), "sort_fun_seq"));
  sort_timings.push_back(std::make_tuple(run_sorting_timer(containers, sort_fun_seq_not_rev), "sort_fun_seq_not_rev"));
  sort_timings.push_back(std::make_tuple(run_sorting_timer(containers, sort_fun_bubblesort), "sort_fun_bubblesort"));

  for (auto [time, sort_algorithm] : sort_timings){
    std::cout << sort_algorithm << " took " << std::get<0>(time) / 1000 << " [microsec] on average, avg kmer -> " << std::get<1>(time) << " [nanosec]" << std::endl; 
  }
}

using out_t = float;

struct new_RI_Kmer {
  Eigen::Array4f next_char_prob;
  // std::array<out_t, 4> next_char_prob;
  int integer_rep;
  bool is_null = true;

  new_RI_Kmer() = default;

  new_RI_Kmer(const vlmc::VLMCKmer &old_kmer){
      Eigen::Array4f tmp = {old_kmer.next_symbol_counts[0], 
                              old_kmer.next_symbol_counts[1],
                              old_kmer.next_symbol_counts[2],
                              old_kmer.next_symbol_counts[3]};
      out_t child_count = tmp.sum() + 4;
      this->next_char_prob = (tmp + pseudo_count_amount) / child_count;
      // out_t child_count = old_kmer.next_symbol_counts[0] + old_kmer.next_symbol_counts[1] + old_kmer.next_symbol_counts[2] + old_kmer.next_symbol_counts[3] + 4;
      // this->next_char_prob = {(old_kmer.next_symbol_counts[0] + pseudo_count_amount) / child_count,
      //                         (old_kmer.next_symbol_counts[1] + pseudo_count_amount) / child_count,
      //                         (old_kmer.next_symbol_counts[2] + pseudo_count_amount) / child_count,
      //                         (old_kmer.next_symbol_counts[3] + pseudo_count_amount) / child_count};
      this->integer_rep = get_index_rep(old_kmer);
      this->is_null = false;
  }

  ~new_RI_Kmer() = default;
    
  new_RI_Kmer(const int temp) : integer_rep{temp} {}

  int get_index_rep(const vlmc::VLMCKmer &kmer) {
    int integer_value = 0;
    int offset = 1;
    for (int i = kmer.length - 1; i >= 0; i--) {
      auto kmer_2_bits = extract2bits(kmer, i) + 1;
      integer_value += (kmer_2_bits * offset);
      offset *= 4;
    }
    return integer_value;
  }
    
  inline uchar extract2bits(const vlmc::VLMCKmer &kmer, uint32 pos) const {
    uchar row = pos >> 5;
    uchar pos_in_row = pos & 31;
    uchar n_shift_pos_to_end = (62 - pos_in_row * 2);
    return (kmer.kmer_data[row] >> n_shift_pos_to_end) & 3;
  }

  inline bool operator<(const new_RI_Kmer &kmer) const {
    return this->integer_rep < kmer.integer_rep;
  };
  inline bool operator>(const new_RI_Kmer &kmer) const {
    return this->integer_rep > kmer.integer_rep;
  };
  inline bool operator>=(const new_RI_Kmer &kmer) const {
    return this->integer_rep >= kmer.integer_rep;
  };
  inline bool operator<=(const new_RI_Kmer &kmer) const {
    return this->integer_rep <= kmer.integer_rep;
  };
  inline bool operator==(const new_RI_Kmer &kmer) const {
    return this->integer_rep == kmer.integer_rep;
  };
};


using tmp_Kmer = new_RI_Kmer;
using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;
#include "utils.hpp"

out_t normalise_dvstar(out_t dot_product, out_t left_norm,
                        out_t right_norm) {

  left_norm = std::sqrt(left_norm);
  right_norm = std::sqrt(right_norm);
  if (left_norm == 0 || right_norm == 0) {
    return 1.0;
  } else {
    out_t Dvstar = dot_product / (left_norm * right_norm);

    out_t dvstar = 0.5 * (1 - Dvstar);

    out_t angular_distance = 2 * std::acos(Dvstar) / M_PI;
    if (isnan(angular_distance)) {
      return 0.0;
    } else {
      return angular_distance;
    }
  }
}

void get_kmer_vector_test(std::filesystem::path path, std::vector<std::vector<tmp_Kmer>>& vlmcs, int max_size){
  std::ifstream ifs(path, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);
  std::vector<tmp_Kmer> vlmc{};
  vlmc::VLMCKmer input_kmer{};
  int count = 0; 

  while (ifs.peek() != EOF){
    if (count > max_size) break; 
    archive(input_kmer);
    tmp_Kmer ri_kmer(input_kmer); 
    vlmc.push_back(ri_kmer);
    count++; 
  }
  ifs.close();
  std::sort(vlmc.begin(), vlmc.end());
  vlmcs.push_back(vlmc);
}

matrix_t iterate_kmers_bench(int max_size){
  std::filesystem::path path_fst{"./data/benchmarking/human/large"};
  std::vector<std::vector<tmp_Kmer>> vlmcs{};
  for (const auto& dir_entry : recursive_directory_iterator(path_fst)) {
    get_kmer_vector_test(dir_entry.path(), vlmcs, max_size);
  }
  
  std::cout << sizeof(vlmcs[0]) + sizeof(tmp_Kmer) * vlmcs[0].capacity() << std::endl;

  matrix_t distances{vlmcs.size(), vlmcs.size()};
  // auto start = std::chrono::steady_clock::now();

  auto fun = [&](auto x, auto y) {
    std::vector<tmp_Kmer> &left_vlmc = vlmcs[x];
    std::vector<tmp_Kmer> &right_vlmc = vlmcs[y];

    out_t dot_product = 0.0;
    out_t left_norm = 0.0;
    out_t right_norm = 0.0; 

    // auto dvstar_fun = [&](auto &left_v, auto &right_v) {
    //   dot_product += (left_v.next_char_prob[0] * right_v.next_char_prob[0]) + (left_v.next_char_prob[1] * right_v.next_char_prob[1]) + (left_v.next_char_prob[2] * right_v.next_char_prob[2]) + (left_v.next_char_prob[3] * right_v.next_char_prob[3]);
    //   left_norm += std::pow(left_v.next_char_prob[0], 2.0) + std::pow(left_v.next_char_prob[1], 2.0) + std::pow(left_v.next_char_prob[2], 2.0) + std::pow(left_v.next_char_prob[3], 2.0);
    //   right_norm += std::pow(right_v.next_char_prob[0], 2.0) + std::pow(right_v.next_char_prob[1], 2.0) + std::pow(right_v.next_char_prob[2], 2.0) + std::pow(right_v.next_char_prob[3], 2.0);
    // };

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

    distances(x,y) = normalise_dvstar(dot_product, left_norm, right_norm); 
  };

  // for (int x = 0; x < vlmcs.size(); x++){
  //   for (int y = 0; y < vlmcs.size(); y++){
  //     fun(x,y); 
  //   }
  // }

  utils::matrix_recursion(0, vlmcs.size(), 0, vlmcs.size(), fun); 

  // std::cout << "Size of RI_kmer struct " << sizeof(RI_Kmer) << " actual object " << sizeof(vlmcs[0][0]) << std::endl; 
  // auto end = std::chrono::steady_clock::now();
  // auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  // std::cout << "Time : " << time << " [micro sec]" << std::endl; 
  // std::cout << "Dot products for sanity check: " << distances(0,0) << "\n";
  // std::cout << "Size of left: " << sizeof(left) + sizeof(RI_Kmer) * left.capacity() << " bytes" << "\n";
  // std::cout << "Size of right: " << sizeof(right) + sizeof(RI_Kmer) * right.capacity() << " bytes" << "\n";
  return distances;
}

int main(int argc, char *argv[]){
  // int num_items = 1500;

  // run_timer<container::VLMC_vector>("Vector");
  //run_timer<container::VLMC_Indexing>("Indexing");
  //run_timer<container::VLMC_sorted_vector>("Sorted Vector");
  //run_timer<container::VLMC_B_tree>("B-tree");
  //run_timer<container::VLMC_hashmap>("Hashmap");
  //run_timer<container::VLMC_Combo>("Combo");
  //run_timer<container::VLMC_Veb>("Veb-tree");
  //benchmark_read_in_kmer();
  //benchmark_kmer_comparison();
  //benchmark_container_inv_sqrt();
  //benchmark_kmer_major_load_calc(8);
  //benchmark_calculate_distance_major();
  //print_kmers_to_file();
  int max_size = std::stoi(argv[1]);
  auto array = iterate_kmers_bench(max_size);
  for (int x = 0; x < array.rows(); x++){
    for (int y = 0; y < array.cols(); y++){
      std::cout << array(x,y) << " ";
    }
    std::cout << std::endl; 
  } 
}