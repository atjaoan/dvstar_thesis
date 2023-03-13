#pragma once 

#include <math.h>
#include <map>
#include <unordered_map>

#include "vlmc_container.hpp"
#include "cluster_container.hpp"
#include "read_in_kmer.hpp"
#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_from_kmers/estimators.hpp"
#include "utils.hpp"

namespace distance {

using vlmc_c = container::VLMC_Container;  
using Kmer = container::RI_Kmer;
using bucket_t = std::unordered_multimap<int, container::Kmer_Pair>::local_iterator;
using matrix_t = Eigen::MatrixXd;

std::array<std::array<double, 4>, 2>
get_components(const Kmer &left, const Kmer &left_background,
               const Kmer &right, const Kmer &right_background) {
  auto left_probs = left.next_char_prob;
  auto left_probs_background = left_background.next_char_prob;
  auto right_probs = right.next_char_prob;
  auto right_probs_background = right_background.next_char_prob;

  return {std::array<double, 4>{
              left_probs[0] / std::sqrt(left_probs_background[0]),
              left_probs[1] / std::sqrt(left_probs_background[1]),
              left_probs[2] / std::sqrt(left_probs_background[2]),
              left_probs[3] / std::sqrt(left_probs_background[3])},
          std::array<double, 4>{
              right_probs[0] / std::sqrt(right_probs_background[0]),
              right_probs[1] / std::sqrt(right_probs_background[1]),
              right_probs[2] / std::sqrt(right_probs_background[2]),
              right_probs[3] / std::sqrt(right_probs_background[3])}};
}

std::array<std::array<double, 4>, 2>
get_components(const Kmer &left, const Kmer &right) {
  auto left_probs = left.next_char_prob;
  auto right_probs = right.next_char_prob;

  return {std::array<double, 4>{
              left_probs[0],
              left_probs[1],
              left_probs[2],
              left_probs[3]},
          std::array<double, 4>{
              right_probs[0],
              right_probs[1],
              right_probs[2],
              right_probs[3]}};
}

std::string get_background_context(const std::string &state,
                                   const size_t background_order) {
  if (state.size() <= background_order) {
    return state; // <- This will never happen 
  } else {
    size_t background = state.size() - background_order;
    return state.substr(background);
  }
}

double normalise_dvstar(double dot_product, double left_norm,
                        double right_norm) {
  left_norm = std::sqrt(left_norm);
  right_norm = std::sqrt(right_norm);

  if (left_norm == 0 || right_norm == 0) {
    return 1.0;
  } else {
    double Dvstar = dot_product / (left_norm * right_norm);

    double dvstar = 0.5 * (1 - Dvstar);

    double angular_distance = 2 * std::acos(Dvstar) / M_PI;
    if (isnan(angular_distance)) {
      return 0.0;
    } else {
      return angular_distance;
    }
  }
}

double dvstar(vlmc_c &left, vlmc_c &right, size_t background_order){

  double dot_product = 0.0;

  double left_norm = 0.0;
  double right_norm = 0.0;

  auto dvstar_fun = [&](auto &left_v, auto &right_v) {
    dot_product += (left_v.next_char_prob * right_v.next_char_prob).sum();
    left_norm += left_v.next_char_prob.square().sum();
    right_norm += right_v.next_char_prob.square().sum();
    };

  if (left.size() < right.size()){
    left.iterate_kmers(left, right, dvstar_fun);
  } else {
    right.iterate_kmers(right, left, dvstar_fun);
  }
      
  return normalise_dvstar(dot_product, left_norm, right_norm);
}

void dvstar_kmer_major(bucket_t &left_begin, bucket_t &left_end, bucket_t &right_begin, bucket_t &right_end, 
                      matrix_t &dot_prod, matrix_t &left_norm, matrix_t &right_norm){
  auto size_left = std::distance(left_begin, left_end);
  auto size_right = std::distance(right_begin, right_end);
  auto rec_fun = [&](size_t &left, size_t &right) { 
    auto left_it = std::next(left_begin, left);
    auto right_it = std::next(right_begin, right);
    if (left_it->second.kmer.integer_rep == right_it->second.kmer.integer_rep){
      auto left_id = left_it->second.id;
      auto right_id = right_it->second.id; 
      dot_prod(left_id, right_id) += (left_it->second.kmer.next_char_prob * right_it->second.kmer.next_char_prob).sum();
      left_norm(left_id, right_id) += left_it->second.kmer.next_char_prob.square().sum();
      right_norm(left_id, right_id) += right_it->second.kmer.next_char_prob.square().sum();
    }
  };

  utils::matrix_recursion(0, size_left, 0, size_right, rec_fun);
}
}