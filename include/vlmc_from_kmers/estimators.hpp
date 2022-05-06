#pragma once

#include <cmath>
#include <numeric>

#include "kmer.hpp"

namespace vlmc {

using estimator_f =
    std::function<std::tuple<bool, double>(const VLMCKmer &, const VLMCKmer &)>;

std::array<double, 4>
get_next_symbol_probabilities(const VLMCKmer &node,
                              const double pseudo_count_amount) {
  // pseudo-counts => +4
  double sum =
      std::accumulate(node.next_symbol_counts.begin(),
                      node.next_symbol_counts.end(), pseudo_count_amount * 4);

  return {(double(node.next_symbol_counts[0]) + pseudo_count_amount) / sum,
          (double(node.next_symbol_counts[1]) + pseudo_count_amount) / sum,
          (double(node.next_symbol_counts[2]) + pseudo_count_amount) / sum,
          (double(node.next_symbol_counts[3]) + pseudo_count_amount) / sum};
}

double kl_divergence(const VLMCKmer &child, const VLMCKmer &parent,
                     const double pseudo_count_amount) {
  auto child_probs = get_next_symbol_probabilities(child, pseudo_count_amount);
  auto parent_probs =
      get_next_symbol_probabilities(parent, pseudo_count_amount);

  double value = 0.0;

  for (int i = 0; i < 4; i++) {
    double child_prob = child_probs[i];
    double parent_prob = parent_probs[i];
    value += child_prob * std::log(child_prob / parent_prob);
  }

  value *= child.count;

  return value;
}

double peres_shields_delta(const VLMCKmer &child, const VLMCKmer &parent,
                           const double pseudo_count_amount) {
  auto child_probs = get_next_symbol_probabilities(child, pseudo_count_amount);
  auto parent_probs =
      get_next_symbol_probabilities(parent, pseudo_count_amount);

  double fluctuation = 0.0;

  for (int i = 0; i < 4; i++) {
    double child_prob = child_probs[i];
    double parent_prob = parent_probs[i];

    double delta_i = parent.count * parent_prob - child_prob * parent.count;

    fluctuation = std::max(fluctuation, delta_i);
  }

  return fluctuation;
}

std::function<std::tuple<bool, double>(const VLMCKmer &, const VLMCKmer &)>
peres_shield_estimator(double sequence_length,
                       const double pseudo_count_amount) {
  double threshold = std::pow(sequence_length, 3.0 / 4.0);

  auto fun = [&, threshold, pseudo_count_amount](
                 const VLMCKmer &child,
                 const VLMCKmer &parent) -> std::tuple<bool, double> {
    auto fluctuation = peres_shields_delta(child, parent, pseudo_count_amount);

    return {fluctuation < threshold, fluctuation};
  };
  return fun;
}

std::function<std::tuple<bool, double>(const VLMCKmer &, const VLMCKmer &)>
kl_estimator(const double threshold, const double pseudo_count_amount) {
  auto fun = [&, threshold, pseudo_count_amount](
                 const VLMCKmer &child,
                 const VLMCKmer &parent) -> std::tuple<bool, double> {
    auto divergence = kl_divergence(child, parent, pseudo_count_amount);

    return {divergence <= threshold, divergence};
  };
  return fun;
}

} // namespace vlmc
