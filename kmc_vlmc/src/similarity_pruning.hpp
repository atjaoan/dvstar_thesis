#pragma once

#include <array>
#include <filesystem>
#include <numeric>
#include <vector>

#include "kmer.hpp"
#include "kmer_container.hpp"
#include "kmers_per_level.hpp"

// Need to keep track of which kmers have children, as those with children
// can't be removed.
struct PstKmer {
  VLMCKmer kmer;
  bool hasChildren = false;
  bool real_child = false;
};

std::array<double, 4> get_next_symbol_probabilities(VLMCKmer &node) {
  // pseudo-counts => +4
  double sum = std::accumulate(node.next_symbol_counts.begin(),
                               node.next_symbol_counts.end(), 4.0);

  return {double(node.next_symbol_counts[0] + 1) / sum,
          double(node.next_symbol_counts[1] + 1) / sum,
          double(node.next_symbol_counts[2] + 1) / sum,
          double(node.next_symbol_counts[3] + 1) / sum};
}

double kl_divergence(VLMCKmer &child, VLMCKmer &parent) {
  auto child_probs = get_next_symbol_probabilities(child);
  auto parent_probs = get_next_symbol_probabilities(parent);

  double value = 0.0;

  for (int i = 0; i < 4; i++) {
    double child_prob = child_probs[i];
    double parent_prob = parent_probs[i];
    value += child_prob * std::log(child_prob / parent_prob);
  }

  value *= child.count;

  return value;
}

bool process_parent(VLMCKmer &prev_kmer, VLMCKmer &kmer,
                    KMersPerLevel<PstKmer> &kmers_per_level,
                    cereal::BinaryOutputArchive &oarchive,
                    const std::function<bool(double)> &remove_node) {
  auto &children = kmers_per_level[prev_kmer.length];

  int removed_children = 0;
  int n_real_children = 0;
  for (auto &[child, has_children, real_child] : children) {
    if (!real_child) {
      continue;
    }
    n_real_children++;

    if (has_children) {
      // Can't remove nodes with children.
      child.divergence = -1.0;
      //      child.output(std::cout);
      oarchive(child);
      continue;
    }

    auto divergence = kl_divergence(child, kmer);

    if (remove_node(divergence)) {
      removed_children++;
      //      child.divergence = divergence;
      //      child.output(stream);
    } else {
      child.divergence = divergence;
      //      child.output(std::cout);
      oarchive(child);
    }
  }

  // Reset children
  kmers_per_level.reset_depth(kmer.length + 1);

  return removed_children != n_real_children;
}

void similarity_prune(VLMCKmer &prev_kmer, VLMCKmer &kmer,
                      KMersPerLevel<PstKmer> &kmers_per_level,
                      cereal::BinaryOutputArchive &oarchive,
                      const std::function<bool(double)> &remove_node) {
  bool has_children = true;

  if (prev_kmer.length > kmer.length) {
    // Has different length, kmer must be parent of the previous nodes
    has_children =
        process_parent(prev_kmer, kmer, kmers_per_level, oarchive, remove_node);
  } else {
    // Has same length, or shorter so the kmers have to be children
    // of the same node, or the first kmer in a new set of children.
    has_children = false;
  }

  if (kmer.length == 0) {
    // Root node, should be the last node, and always outputted.
    oarchive(kmer);
  } else {
    auto start_char_pos = kmer.char_pos(0);
    kmers_per_level[kmer.length][start_char_pos] = {kmer, has_children, true};
  }
}

template <int kmer_size>
void similarity_pruning(KmerContainer<> &container,
                        cereal::BinaryOutputArchive &oarchive,
                        const std::function<bool(double)> &remove_node) {

  KMersPerLevel<PstKmer> kmers_per_level{4, kmer_size + 1};

  VLMCTranslator kmer_api(kmer_size);
  VLMCKmer prev_kmer = kmer_api.construct_vlmc_kmer();

  // Example container output:
  // TTTTTTTT
  // TTTTTTTG
  // TTTTTTTC
  // TTTTTTTA
  // TTTTTTT
  // TTTTTTGT

  container.for_each([&](VLMCKmer &kmer) {
    similarity_prune(prev_kmer, kmer, kmers_per_level, oarchive, remove_node);

    prev_kmer = kmer;
  });
}
