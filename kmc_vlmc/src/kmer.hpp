#pragma once

#include <kmc_file.h>

#include <array>
#include <bitset>

#include <cereal/cereal.hpp>
#include <cereal/types/array.hpp>

// For stxxl, this class needs to be POD
// Otherwise, borrows implementations from KMCs CKmerApi,
// but disregards the byte alignment.
struct VLMCKmer {
  VLMCKmer() = default;
  VLMCKmer(uint32 length_, size_t count_,
           std::array<size_t, 4> next_symbol_counts_)
      : length(length_), count(count_), next_symbol_counts(next_symbol_counts_),
        divergence(-1.0) {
    this->n_rows = (length_ / 32) + 1;
  }

  ~VLMCKmer() = default;

  std::array<uint64, 2> kmer_data;
  size_t count;
  std::array<size_t, 4> next_symbol_counts;
  double divergence;
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

  // This method lets cereal know which data members to serialize
  template <class Archive> void serialize(Archive &archive) {
    archive(kmer_data, length, n_rows, count, next_symbol_counts, divergence);
  }

  /**
   * Finds the first position between the two kmers that is different.
   * E.g. for two kmers "ACGT" and "ACTT", it would return 2
   *
   * @param current_kmer
   * @param prev_kmer
   * @return position the k-mers differ at.
   */
  static int get_first_differing_position(VLMCKmer &current_kmer,
                                          VLMCKmer &prev_kmer,
                                          const int current_kmer_prefix = 0,
                                          const int prev_kmer_prefix = 0) {
    auto n_rows = current_kmer.n_rows;
    int offset = 0;

    for (uint32 row_counter = 0; row_counter < n_rows; row_counter++) {
      unsigned long long current_data = current_kmer.kmer_data[row_counter];
      unsigned long long prev_data = prev_kmer.kmer_data[row_counter];

      if (row_counter == 0) {
        current_data = current_data << (current_kmer_prefix * 2);
        prev_data = prev_data << (prev_kmer_prefix * 2);
      }

      unsigned long long diff = (current_data) ^ (prev_data);

      // GCC instruction, counts leading zeros.
      int n_leading_zeros = __builtin_clzll(diff);

      if (n_leading_zeros != 63) {
        int diff_pos = n_leading_zeros / 2;
        int final_diff_pos = diff_pos + offset;

        if (final_diff_pos >= current_kmer.length - current_kmer_prefix) {
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
    // pos >> 5 = pos // 31
    // pos & 31 = remainder(pos, 31)

    int row = pos >> 5;
    int pos_in_row = pos & 31;
    int n_shift_pos_to_end = (62 - pos_in_row * 2);
    return (kmer_data[row] >> n_shift_pos_to_end) & 3;
  }

  inline int char_pos(int pos) const {
    auto bits = extract2bits(pos);

    return bits;
  }

  inline bool reverse_less_than(const VLMCKmer &kmer) const {
    int this_pos = this->length - 1;
    int kmer_pos = kmer.length - 1;

    while (this_pos >= 0 && kmer_pos >= 0) {
      auto this_2_bits = this->extract2bits(this_pos--);
      auto kmer_2_bits = kmer.extract2bits(kmer_pos--);

      if (this_2_bits != kmer_2_bits) {
        return this_2_bits > kmer_2_bits;
      }
    }

    if (kmer.length != this->length) {
      return kmer.length < this->length;
    }
    return kmer.count > this->count;
  }

  inline bool operator<(const VLMCKmer &kmer) const {
    int min_length = std::min(kmer.length, length);

    for (int i = 0; i < min_length; i++) {
      auto this_2_bits = this->extract2bits(i);
      auto kmer_2_bits = kmer.extract2bits(i);
      if (this_2_bits != kmer_2_bits) {
        return this_2_bits < kmer_2_bits;
      }
    }
    return kmer.length > this->length;
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
    stream << divergence << std::endl;
  }
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
    } else if (this->kmer_length > 0) {
      new_kmer.kmer_data[0] = this->kmer_data[0] << (this->byte_alignment) * 2;
    } else {
      new_kmer.kmer_data[0] = 0;
    }
    return new_kmer;
  }
};

template <int MAX_K> struct VirtualKMerComparator {
  virtual bool operator()(const VLMCKmer &a, const VLMCKmer &b) const {
    return false;
  };
  virtual VLMCKmer min_value() const { return VLMCKmer(0, 0, {}); }
  virtual VLMCKmer max_value() const { return VLMCKmer(0, 0, {}); }
};

template <int MAX_K>
struct KMerComparator : public VirtualKMerComparator<MAX_K> {
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

template <int MAX_K>
struct ReverseKMerComparator : public VirtualKMerComparator<MAX_K> {
  bool operator()(const VLMCKmer &a, const VLMCKmer &b) const {
    return a.reverse_less_than(b);
  }
  VLMCKmer min_value() const {
    VLMCKmer max_kmer(MAX_K, (size_t)-1, {});
    for (int row = 0; row < max_kmer.n_rows; row++) {
      max_kmer.kmer_data[row] = 0xFFFFFFFFFFFFFFFF;
    }
    return max_kmer;
  }
  VLMCKmer max_value() const {
    std::array<size_t, 4> vec{};
    return VLMCKmer(0, 0, vec);
  }
};

template <int MAX_K> struct KMerReverseKeyExtractor {
  typedef uint64 key_type;
  uint64 upper_bound = (1 - std::pow(4, MAX_K + 1)) / (1 - 4);

  uint64 n_kmers_with_length(uint32_t len) const {
    return (1 - std::pow(4, len + 1)) / (1 - 4);
  }

  key_type operator()(const VLMCKmer &a) const {
    uint64 key = 0;

    for (int i = 0; i < a.length; i++) {
      int bits = a.extract2bits(a.length - i - 1);

      // Length is reversed - suffixes should count as having only the most
      // significant digits.
      int length_from_end = MAX_K - i - 1;

      // For every position we add the number of k-mers with a specific
      // length that have to come before this one, multiplied by the bit values.
      // This can be viewed as dividing the space of all MAX_K k-mers
      // into divisions for every suffix, where every such division at a
      // specific length is n_kmers_with_length(i) large.  Thus, we need
      // to find which such division this kmer belongs to.
      // Plus one so e.g. AT and T don't get same values.

      key += bits * this->n_kmers_with_length(length_from_end) + 1;
    }

    return upper_bound - key;
  }

  VLMCKmer min_value() const {
    VLMCKmer max_kmer(MAX_K, (size_t)-1, {});
    for (int row = 0; row < max_kmer.n_rows; row++) {
      max_kmer.kmer_data[row] = 0xFFFFFFFFFFFFFFFF;
    }
    return max_kmer;
  }
  VLMCKmer max_value() const {
    std::array<size_t, 4> vec{};
    return VLMCKmer(0, 0, vec);
  }
};
