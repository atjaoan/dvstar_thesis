#pragma once 

//temp
#include "read_in_kmer.hpp"

#include "vlmc_container.hpp"
#include "read_in_kmer.hpp"
#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_from_kmers/estimators.hpp"

namespace distance {

using vlmc_c = container::VLMC_Container;  
using Kmer = container::RI_Kmer;

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

std::string get_background_context(const std::string &state,
                                   const size_t background_order) {
  if (state.size() <= background_order) {
    return state;
  } else {
    size_t background = state.size() - background_order;
    return state.substr(background);
  }
}

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
    auto right_kmer = right_kmers.find(left_kmer);
    if (right_kmer.is_null){
      f_not_shared(left_kmer, right_kmer);
    } else {
      f(left_kmer, right_kmer);
    }
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
    if (std::isnan(angular_distance)) {
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

  iterate_kmers(
      left, right, [&](const Kmer &left_v, const Kmer &right_v) {
        if (left_v.length <= background_order) {
          return;
        }
        // const auto background_context = get_background_context(left_v.to_string(), background_order);
        std::string background_context{""};

        //Create new kmer with new context, should be reconsidered
        vlmc::VLMCTranslator kmer{static_cast<int>(background_context.size())};
        if (!background_context.empty()) {
          kmer.from_string(background_context);
        }
        Kmer context = kmer.construct_vlmc_kmer();
        container::RI_Kmer ri_context{context};

        auto left_kmer_background = left.find(ri_context);
        auto right_kmer_background = right.find(ri_context);

        auto [left_comp, right_comp] = get_components(
            left_v, left_kmer_background, right_v, right_kmer_background);

        for (int i = 0; i < 4; i++) { 
          dot_product += left_comp[i] * right_comp[i];
          left_norm += std::pow(left_comp[i], 2.0);
          right_norm += std::pow(right_comp[i], 2.0);
        }
      },
      [](auto l, auto r) {});

  return normalise_dvstar(dot_product, left_norm, right_norm);
}
}