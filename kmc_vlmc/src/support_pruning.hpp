#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <string>

#include <kmc_file.h>

int get_first_differing_position(CKmerAPI &current_kmer, CKmerAPI &prev_kmer) {
  auto no_of_rows = current_kmer.no_of_rows;
  int offset = - current_kmer.byte_alignment;

  for (uint32 row_counter = 0; row_counter < no_of_rows; row_counter++) {
    unsigned long long current_data = current_kmer.kmer_data[row_counter];
    unsigned long long prev_data = prev_kmer.kmer_data[row_counter];

    unsigned long long diff = (current_data) ^ (prev_data);

    // GCC instruction, counts leading zeros.
    int n_leading_zeros = __builtin_clzll(diff);

    if (n_leading_zeros != 63) {
      int diff_pos = n_leading_zeros / 2;
      int final_diff_pos = diff_pos + offset;

      if (final_diff_pos > current_kmer.kmer_length) {
        return -1;
      } else {
        return final_diff_pos;
      }
    }
    offset += 32;
  }
  return -1;
}

void process_kmer(CKmerAPI &current_kmer, CKmerAPI &prev_kmer) {
  std::cout << current_kmer.to_string() << " ";
  std::cout << get_first_differing_position(current_kmer, prev_kmer);

  std::cout << std::endl;
}

