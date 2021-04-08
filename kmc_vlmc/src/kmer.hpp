#pragma once

#include <kmc_file.h>

#include <array>
#include <bitset>

#include <stxxl/sorter>

// For stxxl, this class needs to be POD
// Otherwise, borrows implementations from KMCs CKmerApi,
// but disregards the byte alignment.
class VLMCKmer {
public:
  VLMCKmer() = default;
  VLMCKmer(uint32 length_, size_t count_,
           std::array<size_t, 4> next_symbol_counts_)
      : length(length_), count(count_),
        next_symbol_counts(next_symbol_counts_) {
    this->n_rows = (length_ / 32) + 1;
  }

  ~VLMCKmer() = default;

  std::array<uint64, 2> kmer_data;
  size_t count;
  std::array<size_t, 4> next_symbol_counts;
  int length;
  int n_rows;

  static constexpr char char_codes[4] = {'A', 'C', 'G', 'T'};

  static VLMCKmer create_prefix_kmer(VLMCKmer &kmer, uint32 length,
                                     size_t count,
                                     std::array<size_t, 4> child_counts) {
    VLMCKmer prefix_kmer{length, count, child_counts};

    for (int row_counter = 0; row_counter < prefix_kmer.n_rows; row_counter++) {
      unsigned long long kmer_data = kmer.kmer_data[row_counter];

      int positions_to_remove_from_row = (32 - length) * 2 % 32;
      if (length > 32 * (row_counter + 1)) {
        positions_to_remove_from_row = 0;
      }
      unsigned long long mask = 0xFFFFFFFFFFFFFFFF
                                << positions_to_remove_from_row;

      prefix_kmer.kmer_data[row_counter] = (kmer_data & mask);
    }

    return prefix_kmer;
  }

  static int get_first_differing_position(VLMCKmer &current_kmer,
                                          VLMCKmer &prev_kmer) {
    auto n_rows = current_kmer.n_rows;
    int offset = 0;

    for (uint32 row_counter = 0; row_counter < n_rows; row_counter++) {
      unsigned long long current_data = current_kmer.kmer_data[row_counter];
      unsigned long long prev_data = prev_kmer.kmer_data[row_counter];

      unsigned long long diff = (current_data) ^ (prev_data);

      // GCC instruction, counts leading zeros.
      int n_leading_zeros = __builtin_clzll(diff);

      if (n_leading_zeros != 63) {
        int diff_pos = n_leading_zeros / 2;
        int final_diff_pos = diff_pos + offset;

        if (final_diff_pos > current_kmer.length) {
          return -1;
        } else {
          return final_diff_pos;
        }
      }
      offset += 32;
    }
    return -1;
  }

  inline uchar extract2bits(uint32 pos) const {
    return (kmer_data[pos >> 5] >> (62 - ((pos & 31) * 2))) & 3;
  }

  inline int char_pos(int pos) const {
    auto bits = extract2bits(pos);

    return bits;
  }

  inline bool operator<(const VLMCKmer &kmer) const {
    int min_length = std::min(kmer.length, length);

    for (int i = 0; i < min_length; i++) {
      auto local_2_bits = this->extract2bits(i);
      auto kmer_2_bits = kmer.extract2bits(i);
      if (local_2_bits != kmer_2_bits) {
        return local_2_bits < kmer_2_bits;
      }
    }
    return kmer.length > length;
  };

  inline std::string to_string() {
    std::string out_string(this->length, ' ');

    uchar *byte_ptr;
    uchar c;
    uint32 cur_string_size = 0;

    for (uint32 row_counter = 0; row_counter < this->n_rows; row_counter++) {
      byte_ptr = reinterpret_cast<uchar *>(&kmer_data[row_counter]);

      byte_ptr += 7; // Read from the left

      for (uint32 i = 0; (i < length) && (i < 32); i += 4) {

        c = 0xc0 & *byte_ptr; // 11000000
        c = c >> 6;
        out_string[cur_string_size++] = char_codes[c];

        if (cur_string_size == length)
          break;

        c = 0x30 & *byte_ptr; // 00110000
        c = c >> 4;
        out_string[cur_string_size++] = char_codes[c];
        if (cur_string_size == length)
          break;

        c = 0x0c & *byte_ptr; // 00001100
        c = c >> 2;
        out_string[cur_string_size++] = char_codes[c];
        if (cur_string_size == length)
          break;

        c = 0x03 & *byte_ptr; // 00000011
        out_string[cur_string_size++] = char_codes[c];
        if (cur_string_size == length)
          break;

        byte_ptr--;
      }
    }

    return out_string;
  }

  void output(std::ostream &stream) {
    stream << this->to_string() << " " << this->count << " ";
    for (auto &c : this->next_symbol_counts) {
      stream << c << " ";
    }
    stream << std::endl;
  }


protected:
};

class VLMCTranslator : public CKmerAPI {
public:
  VLMCTranslator(int length) : CKmerAPI(length) {}
  ~VLMCTranslator() = default;

  VLMCKmer construct_vlmc_kmer() {
    // To make sure the VLMCKmer class stays a POD, we're only allowing k-mers
    // up to 64 or so.  Should be fine.
    VLMCKmer new_kmer{this->kmer_length, 0, {}};

    new_kmer.kmer_data[0] = 0;
    new_kmer.kmer_data[1] = 0;

    if (this->kmer_length > 32 - this->byte_alignment) {
      for (int pos = 0; pos < this->kmer_length; pos++) {
        uint32 row = pos >> 5;
        uchar val = this->extract2bits(pos);
        new_kmer.kmer_data[row] += (uint64)val << (62 - ((pos & 31) * 2));
      }
    } else if (this->kmer_length > 0){
      new_kmer.kmer_data[0] = this->kmer_data[0] << (this->byte_alignment) * 2;
    }
    return new_kmer;
  }
};

template <int MAX_K> struct KMerComparator {
  bool operator()(const VLMCKmer &a, const VLMCKmer &b) const { return a < b; }
  VLMCKmer min_value() const {
    std::array<size_t, 4> vec{};
    return VLMCKmer(0, 0, vec);
  }
  VLMCKmer max_value() const {
    VLMCKmer max_kmer(MAX_K, 0, {});
    for (int row = 0; row < max_kmer.n_rows; row++) {
      max_kmer.kmer_data[row] = 0xFFFFFFFFFFFFFFFF;
    }
    return max_kmer;
  }
};

template <int MAX_K> struct ReverseKMerComparator {
  bool operator()(const VLMCKmer &a, const VLMCKmer &b) const {
    // TODO this needs to compare the kmers backwards.
    return a < b;
  }
  VLMCKmer min_value() const {
    std::array<size_t, 4> vec{};
    return VLMCKmer(0, 0, vec);
  }
  VLMCKmer max_value() const {
    VLMCKmer max_kmer(MAX_K, 0, {});
    for (int row = 0; row < max_kmer.n_rows; row++) {
      max_kmer.kmer_data[row] = 0xFFFFFFFFFFFFFFFF;
    }
    return max_kmer;
  }
};

template <int MAX_K>
using kmer_sorter =
    stxxl::sorter<VLMCKmer, KMerComparator<MAX_K>, 1 * 1024 * 1024>;
