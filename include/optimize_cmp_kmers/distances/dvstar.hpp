#pragma once 

#include <math.h>
#include <map>
#include <unordered_map>

#include "vlmc_container.hpp"
#include "read_in_kmer.hpp"
#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_from_kmers/estimators.hpp"

namespace distance {

using vlmc_c = container::VLMC_Container;  
using Kmer = container::RI_Kmer;
using kmer_c = container::Kmer_Cluster;

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

/*
// Perhaps remove this since it is now placed in the vlmc_container.hpp script

void iterate_kmers(
    vlmc_c &left_kmers, vlmc_c &right_kmers,
    const std::function<void(const Kmer &left, const Kmer &right)> &f,
    const std::function<void(const Kmer &left, const Kmer &right)>
        &f_not_shared) {
  for (size_t i = left_kmers.get_min_kmer_index() ; i <= left_kmers.get_max_kmer_index(); i++) {
    const Kmer &left_kmer = left_kmers.get(i);
    if (left_kmer.is_null){
      continue; 
    }
    auto right_kmer = right_kmers.find(left_kmer.integer_rep);
    if (right_kmer.is_null){
      f_not_shared(left_kmer, right_kmer);
    } else {
      f(left_kmer, right_kmer);
    }
  }
}
*/

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
    //for (int i = 0; i < 4; i++) { 
    //  dot_product += left_v.next_char_prob[i] * right_v.next_char_prob[i];
    //  left_norm += std::pow(left_v.next_char_prob[i], 2.0);
    //  right_norm += std::pow(right_v.next_char_prob[i], 2.0);
    //}
    };

  if (left.size() < right.size()){
    left.iterate_kmers(left, right, dvstar_fun);
  } else {
    right.iterate_kmers(right, left, dvstar_fun);
  }
      
  return normalise_dvstar(dot_product, left_norm, right_norm);
}

double dvstar(bucket_t &left, bucket_t &right, size_t background_order){

  double dot_product = 0.0;

  double left_norm = 0.0;
  double right_norm = 0.0;

  auto dvstar_fun = [&](auto &left_v, auto &right_v) {
    dot_product += (left_v.next_char_prob * right_v.next_char_prob).sum();
    left_norm += left_v.next_char_prob.square().sum();
    right_norm += right_v.next_char_prob.square().sum();
    //for (int i = 0; i < 4; i++) { 
    //  dot_product += left_v.next_char_prob[i] * right_v.next_char_prob[i];
    //  left_norm += std::pow(left_v.next_char_prob[i], 2.0);
    //  right_norm += std::pow(right_v.next_char_prob[i], 2.0);
    //}
    };

  if (left.size() < right.size()){
    left.iterate_kmers(left, right, dvstar_fun);
  } else {
    right.iterate_kmers(right, left, dvstar_fun);
  }
      
  return normalise_dvstar(dot_product, left_norm, right_norm);
}
}