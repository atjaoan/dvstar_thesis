#include <Eigen/Dense>
#include <chrono> 

#include "cluster_container.hpp"
#include "vlmc_container.hpp"
#include "parallel.hpp"
#include "distances/dvstar.hpp"

namespace calculate {

using matrix_t  = Eigen::MatrixXd;
using vlmc_c    = container::VLMC_Container;  
using kmer_pair = container::Kmer_Pair;

template <typename VC>
void calculate_reduced_slice(size_t start_index, size_t stop_index, matrix_t &distances,
                     container::Cluster_Container<VC> &cluster_left, container::Cluster_Container<VC> &cluster_right,
                     const std::function<double(vlmc_c &, vlmc_c &)> &fun) {
  size_t y_bound = 1;
  for (size_t i = start_index; i < stop_index; i++) {
    for (size_t j = y_bound; j < cluster_right.size(); j++) {
      distances(i, j) = fun(cluster_left.get(i), cluster_right.get(j));
    }
    y_bound++;
  }
}

template <typename VC>
void calculate_full_slice_rec(size_t start_index_left, size_t stop_index_left, size_t start_index_right, size_t stop_index_right, matrix_t &distances,
                     container::Cluster_Container<VC> &cluster_left, container::Cluster_Container<VC> &cluster_right,
                     const std::function<double(vlmc_c &, vlmc_c &)> &fun) {
  auto diff_left = stop_index_left - start_index_left;
  auto diff_right = stop_index_right - start_index_right;
  if (diff_left == 1 && diff_right == 1){
    distances(start_index_left, start_index_right) = fun(cluster_left.get(start_index_left), cluster_right.get(start_index_right));
  } else if (diff_right > diff_left){
    auto new_right_index = (stop_index_right + start_index_right) / 2;
    calculate_full_slice_rec(start_index_left, stop_index_left, start_index_right, new_right_index, distances, cluster_left, cluster_right, fun);
    calculate_full_slice_rec(start_index_left, stop_index_left, new_right_index, stop_index_right, distances, cluster_left, cluster_right, fun);
  } else {
    auto new_left_index = (stop_index_left + start_index_left) / 2;
    calculate_full_slice_rec(start_index_left, new_left_index, start_index_right, stop_index_right, distances, cluster_left, cluster_right, fun);
    calculate_full_slice_rec(new_left_index, stop_index_left, start_index_right, stop_index_right, distances, cluster_left, cluster_right, fun);
  }
}

template <typename VC>
void calculate_full_slice(size_t start_index, size_t stop_index, matrix_t &distances,
                     container::Cluster_Container<VC> &cluster_left, container::Cluster_Container<VC> &cluster_right,
                     const std::function<double(vlmc_c &, vlmc_c &)> &fun) {
  calculate_full_slice_rec(start_index, stop_index, 0, cluster_right.size(), distances, cluster_left, cluster_right, fun); 
}

// Inter-directory distances
template <typename VC> 
matrix_t calculate_distances(
    container::Cluster_Container<VC> &cluster, std::function<double(vlmc_c &, vlmc_c &)> &distance_function,
    size_t requested_cores){

  matrix_t distances{cluster.size(), cluster.size()};

  auto fun = [&](size_t start_index, size_t stop_index) {
    calculate_reduced_slice<VC>(start_index, stop_index, distances,
                           cluster, cluster, distance_function);
  };
  parallel::parallelize(cluster.size(), fun, requested_cores);
  return distances; 
}

// For two different dirs
template <typename VC>
matrix_t calculate_distances(
    container::Cluster_Container<VC> &cluster_left, container::Cluster_Container<VC> &cluster_right,
    std::function<double(vlmc_c &, vlmc_c &)> &distance_function,
    size_t requested_cores){

      matrix_t distances{cluster_left.size(), cluster_right.size()};

      auto fun = [&](size_t start_index, size_t stop_index) {
      calculate_full_slice<VC>(start_index, stop_index, std::ref(distances),
                           std::ref(cluster_left), std::ref(cluster_right), distance_function);
      };
      parallel::parallelize(cluster_left.size(), cluster_right.size(), fun, requested_cores);
      // parallel::sequential(cluster_left.size(), cluster_right.size(), fun, requested_cores);
      return distances; 
}

void calculate_kmer_buckets(size_t start_bucket, size_t stop_bucket, 
    matrix_t &distances, matrix_t &dot_prod, matrix_t &left_norm, matrix_t &right_norm,
    container::Kmer_Cluster &cluster_left, container::Kmer_Cluster &cluster_right) {
  for (size_t i = start_bucket; i < stop_bucket; i++) {
    if (cluster_left.is_bucket_empty(i) || cluster_right.is_bucket_empty(i)){
      continue; 
    }
    auto vec_left = cluster_left.get(i);
    auto vec_right = cluster_right.get(i);
    distance::dvstar_kmer_major(vec_left, vec_right, dot_prod, left_norm, right_norm);
  } 
}

void normalize_rec(size_t start_index_left, size_t stop_index_left, size_t start_index_right, size_t stop_index_right, 
  matrix_t &distances, matrix_t &dot_prod, matrix_t &left_norm, matrix_t &right_norm){
  auto diff_left = stop_index_left - start_index_left;
  auto diff_right = stop_index_right - start_index_right;
  if (diff_left == 1 && diff_right == 1){
    distances(start_index_left, start_index_right) = distance::normalise_dvstar(dot_prod(start_index_left, start_index_right), 
                              left_norm(start_index_left, start_index_right), right_norm(start_index_left, start_index_right));
  } else if (diff_right > diff_left){
    auto new_right_index = (stop_index_right + start_index_right) / 2;
    normalize_rec(start_index_left, stop_index_left, start_index_right, new_right_index, distances, dot_prod, left_norm, right_norm);
    normalize_rec(start_index_left, stop_index_left, new_right_index, stop_index_right, distances, dot_prod, left_norm, right_norm);
  } else {
    auto new_left_index = (stop_index_left + start_index_left) / 2;
    normalize_rec(start_index_left, new_left_index, start_index_right, stop_index_right, distances, dot_prod, left_norm, right_norm);
    normalize_rec(new_left_index, stop_index_left, start_index_right, stop_index_right, distances, dot_prod, left_norm, right_norm);
  }
}

matrix_t calculate_distance_major(
    container::Kmer_Cluster &cluster_left, container::Kmer_Cluster &cluster_right,
    size_t requested_cores){

  matrix_t distances = matrix_t::Zero(cluster_left.size(), cluster_right.size());
  matrix_t dot_prod = matrix_t::Zero(cluster_left.size(), cluster_right.size());
  matrix_t left_norm = matrix_t::Zero(cluster_left.size(), cluster_right.size());
  matrix_t right_norm = matrix_t::Zero(cluster_left.size(), cluster_right.size());

  auto fun = [&](size_t start_bucket, size_t stop_bucket) {
    calculate_kmer_buckets(start_bucket, stop_bucket, distances, dot_prod, left_norm, right_norm, cluster_left, cluster_right);
  };
  parallel::parallelize(cluster_left.bucket_count(), fun, requested_cores);


  auto norm_fun = [&](size_t start_vlmc, size_t stop_vlmc) {
    normalize_rec(start_vlmc, stop_vlmc, 0, cluster_right.size(), distances, dot_prod, left_norm, right_norm);
  };

  parallel::parallelize(cluster_left.size(), norm_fun, requested_cores); 

  return distances; 
}

matrix_t calculate_distance_major(
    container::Kmer_Cluster &cluster,
    size_t requested_cores){

  matrix_t distances = matrix_t::Zero(cluster.size(), cluster.size());
  matrix_t dot_prod = matrix_t::Zero(cluster.size(), cluster.size());
  matrix_t left_norm = matrix_t::Zero(cluster.size(), cluster.size());
  matrix_t right_norm = matrix_t::Zero(cluster.size(), cluster.size());

  auto fun = [&](size_t start_bucket, size_t stop_bucket) {
    calculate_kmer_buckets(start_bucket, stop_bucket, distances, dot_prod, left_norm, right_norm, cluster, cluster);
  };
  parallel::parallelize(cluster.bucket_count(), fun, requested_cores);


  auto norm_fun = [&](size_t start_vlmc, size_t stop_vlmc) {
    normalize_rec(start_vlmc, stop_vlmc, 0, cluster.size(), distances, dot_prod, left_norm, right_norm);
  };

  parallel::parallelize(cluster.size(), norm_fun, requested_cores); 

  return distances;
}

}