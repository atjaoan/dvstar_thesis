#pragma once 

#include <functional>
#include <memory>
#include <filesystem>
#include <numeric>

#include "vlmc_from_kmers/kmer.hpp"

/*
  Stores VLMC (multiple k-mers) in a container. 
*/

namespace container{
 
struct RI_Kmer{
    //Could be 3 and infer the fourth probability,
    //Using 4 since it may be useful in SIMD instructions
    std::array<double,4> next_char_prob{};
    int integer_rep;

    RI_Kmer() = default;
    RI_Kmer(const vlmc::VLMCKmer old_kmer){
        this->integer_rep = get_index_rep(old_kmer);
        double child_count = std::accumulate(old_kmer.next_symbol_counts.begin(), old_kmer.next_symbol_counts.end(), 0.0);
        for (size_t i = 0; i < 4; i++){
            this->next_char_prob[i] = old_kmer.next_symbol_counts[i] / child_count;
        }
    }
    ~RI_Kmer() = default;

    int get_index_rep(const vlmc::VLMCKmer &kmer) {
      int integer_value = 0;
      int offset = 1;
      for (int i = kmer.length - 1; i >= 0; i--) {
        auto kmer_2_bits = kmer.extract2bits(i) + 1;
        integer_value += (kmer_2_bits * offset);
        offset *= 4;
      }
      return integer_value;
    }
};
}
