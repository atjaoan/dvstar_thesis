#pragma once 

#include <functional>
#include <memory>
#include <filesystem>
#include <numeric>

#include "vlmc_from_kmers/kmer.hpp"

/*
  Stores VLMC (multiple k-mers) in a container. 
*/

constexpr double pseudo_count_amount = 1.0; 

namespace container{
 
struct RI_Kmer{
    //Could be 3 and infer the fourth probability,
    //Using 4 since it may be useful in SIMD instructions
    size_t length = 0; 
    int integer_rep;
    // int background_rep; <- Should be implemented 
    std::array<double,4> next_char_prob{};
    std::array<uint64, 4> bit_representation;
    bool is_null = true;

    RI_Kmer() = default;
    RI_Kmer(const vlmc::VLMCKmer &old_kmer){
        this->length = old_kmer.length;
        double child_count = std::accumulate(old_kmer.next_symbol_counts.begin(), old_kmer.next_symbol_counts.end(), pseudo_count_amount * 4);
        for (size_t i = 0; i < 4; i++){
            this->next_char_prob[i] = (static_cast<double>(old_kmer.next_symbol_counts[i]) + pseudo_count_amount) / child_count; // Perhaps change static cast to double
            this->bit_representation[i] = old_kmer.kmer_data[i];
        }
        this->integer_rep = get_index_rep(old_kmer);
        this->is_null = false;
        // this->background_rep = background_order_index(this->integer_rep, 0);  <- Should be implemented 
    }
    ~RI_Kmer() = default;

    int get_index_rep(const vlmc::VLMCKmer &kmer) {
      int integer_value = 0;
      int offset = 1;
      for (int i = kmer.length - 1; i >= 0; i--) {
        auto kmer_2_bits = extract2bits(i) + 1;
        integer_value += (kmer_2_bits * offset);
        offset *= 4;
      }
      return integer_value;
    }

    uchar extract2bits(uint32 pos) const {
    // pos >> 5 == pos // 31
    // pos & 31 == remainder(pos, 31)

    uchar row = pos >> 5;
    uchar pos_in_row = pos & 31;
    uchar n_shift_pos_to_end = (62 - pos_in_row * 2);
    return (bit_representation[row] >> n_shift_pos_to_end) & 3;
  }

  int background_order_index(int integer_rep, int order){
    if(integer_rep < std::pow(4, order)) return integer_rep;
    int back_rep = 0;
    int i = 1;
    for(int o = 0; o < order; o++){
      int r = integer_rep % 4;
      if(r == 0) r = 4;
      integer_rep = (integer_rep - r) / 4;
      back_rep += r*i;
      i *= 4;
    }
    return back_rep;
  }

  inline bool operator<(const RI_Kmer &kmer) const {
    return this->integer_rep < kmer.integer_rep;
  };

  inline bool operator==(const RI_Kmer &kmer) const {
    return this->integer_rep == kmer.integer_rep;
  };
};
}
