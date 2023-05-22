#pragma once 

#include <math.h>

#include "vlmc_container.hpp"
#include "cluster_container.hpp"
#include "read_in_kmer.hpp"
#include "utils.hpp"
#include "global_aliases.hpp"

namespace distance {

using RI_Kmer = container::RI_Kmer;
using bucket_t = std::vector<container::Kmer_Pair>;

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

template <typename VC>
out_t dvstar(VC &left, VC &right){
  if (left.size() < right.size()){
    return container::iterate_kmers(left, right);
  } else {
    return container::iterate_kmers(right, left);
  }
}

void dvstar_kmer_major(bucket_t &left_vector, bucket_t &right_vector, 
                      matrix_t &dot_prod, matrix_t &left_norm, matrix_t &right_norm){
  auto rec_fun = [&](size_t &left, size_t &right) { 
    auto left_id = left_vector[left].id;
    auto right_id = right_vector[right].id; 
    dot_prod(left_id, right_id) += (left_vector[left].kmer.next_char_prob * right_vector[right].kmer.next_char_prob).sum();
    left_norm(left_id, right_id) += left_vector[left].kmer.next_char_prob.square().sum();
    right_norm(left_id, right_id) += right_vector[right].kmer.next_char_prob.square().sum();
  };

  utils::matrix_recursion(0, left_vector.size(), 0, right_vector.size(), rec_fun);
}

void dvstar_kmer_major_single(bucket_t &left_vector, bucket_t &right_vector, 
                      matrix_t &dot_prod, matrix_t &left_norm, matrix_t &right_norm){
                        
  auto rec_fun = [&](size_t &left, size_t &right) { 
    if(left_vector[left].id >= right_vector[right].id) return;
    if (left_vector[left].kmer.integer_rep == right_vector[right].kmer.integer_rep){
      auto left_id = left_vector[left].id;
      auto right_id = right_vector[right].id; 
      dot_prod(left_id, right_id) += (left_vector[left].kmer.next_char_prob * right_vector[right].kmer.next_char_prob).sum();
      left_norm(left_id, right_id) += left_vector[left].kmer.next_char_prob.square().sum();
      right_norm(left_id, right_id) += right_vector[right].kmer.next_char_prob.square().sum();
    }
  };

  utils::matrix_recursion(0, left_vector.size(), 0, right_vector.size(), rec_fun);
}
}