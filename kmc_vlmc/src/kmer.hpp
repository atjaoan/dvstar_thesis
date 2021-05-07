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
  VLMCKmer(uint32 length_, uint64 count_,
           std::array<uint64, 4> next_symbol_counts_)
      : length(length_), count(count_), next_symbol_counts(next_symbol_counts_),
        divergence(-1.0), kmer_data() {
    this->n_rows = (length_ / 32.0) + 1; // std::ceil, but no float conversion
  }

  ~VLMCKmer() = default;

  std::array<uint64, 2> kmer_data;
  uint64 count;
  std::array<uint64, 4> next_symbol_counts;
  double divergence;
  uint32 length;
  uint32 n_rows;

  static constexpr char char_codes[4] = {'A', 'C', 'G', 'T'};

  static VLMCKmer create_prefix_kmer(const VLMCKmer &kmer, uint32 length,
                                     uint64 count,
                                     std::array<uint64, 4> child_counts) {
    VLMCKmer prefix_kmer{length, count, child_counts};

    prefix_kmer.kmer_data[0] = 0;
    prefix_kmer.kmer_data[1] = 0;

    uint32 row = length >> 5;

    if (row == 1) {
      prefix_kmer.kmer_data[0] = kmer.kmer_data[0];
    }

    uint32 length_in_row = length & 31;

    uint32 positions_to_remove_from_row = (31 - length_in_row) * 2;
    unsigned long long mask = 0xFFFFFFFFFFFFFFFF
                              << positions_to_remove_from_row;

    prefix_kmer.kmer_data[row] = (kmer.kmer_data[row] & mask);

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
  static int get_first_differing_position(const VLMCKmer &current_kmer,
                                          const VLMCKmer &prev_kmer,
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

  [[nodiscard]] uchar extract2bits(uint32 pos) const {
    // pos >> 5 = pos // 31
    // pos & 31 = remainder(pos, 31)

    int row = pos >> 5;
    int pos_in_row = pos & 31;
    int n_shift_pos_to_end = (62 - pos_in_row * 2);
    return (kmer_data[row] >> n_shift_pos_to_end) & 3;
  }

  [[nodiscard]] int char_pos(int pos) const {
    auto bits = extract2bits(pos);

    return bits;
  }

  [[nodiscard]] static int n_kmers_with_length(uint32_t len) {
    static std::vector<int> vals = {1, 5, 21, 85, 341, 1365};
    //    return vals[len];

    if (len > vals.size()) {
      return (1 - std::pow(4, len + 1)) / (1 - 4);
    } else {
      return vals[len];
    }
  }

  /*
   * Get the value of the first prefix_length characters.
   *
   * Useful for splitting the kmers into subsets based on their first
   * `prefix_length` characters.
   *
   * WARNING: Assumes the prefix is smaller than 32, as this is our only
   * use case.
   *
   * @param prefix_length length of prefix.
   * @return unique value per prefix.
   */
  [[nodiscard]] uint64 get_prefix_index(int prefix_length) const {
    uint64 key = 0;

    int n_shift_pos_to_end = (32 - prefix_length) * 2;
    return kmer_data[0] >> n_shift_pos_to_end;

    for (int i = 0; i < prefix_length; i++) {
      int bits = extract2bits(i);

      // Length is reversed - suffixes should count as having only the most
      // significant digits.
      int length_from_end = prefix_length - i - 1;

      // For every position we add the number of k-mers with a specific
      // length that have to come before this one, multiplied by the bit values.
      // This can be viewed as dividing the space of all prefix_length k-mers
      // into divisions for every suffix, where every such division at a
      // specific length is n_kmers_with_length(i) large.  Thus, we need
      // to find which such division this kmer belongs to.
      // Plus one so e.g. AT and T don't get same values.

      key += bits * VLMCKmer::n_kmers_with_length(length_from_end) + 1;
    }

    return key;
  }

  [[nodiscard]] bool reverse_less_than(const VLMCKmer &kmer) const {
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

    if (kmer.length != this->length) {
      return kmer.length > this->length;
    }
    return kmer.count < this->count;
  };

  [[nodiscard]] std::string to_string() const {
    std::string out_string(this->length, ' ');

    uchar *byte_ptr;
    uchar c;
    uint32 cur_string_size = 0;

    for (uint32 row_counter = 0; row_counter < this->n_rows; row_counter++) {
      byte_ptr = reinterpret_cast<uchar *>(
          &(const_cast<VLMCKmer *>(this)->kmer_data[row_counter]));

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
  explicit VLMCTranslator(int length) : CKmerAPI(length) {}
  ~VLMCTranslator() = default;

  VLMCKmer construct_vlmc_kmer() {
    // To make sure the VLMCKmer class stays a POD, we're only allowing k-mers
    // up to 64 or so.  Should be fine.
    VLMCKmer new_kmer{this->kmer_length, 0, {}};

    new_kmer.kmer_data[0] = 0;
    new_kmer.kmer_data[1] = 0;

    if (this->kmer_length > 31 - this->byte_alignment) {
      for (int pos = 0; pos < this->kmer_length; pos++) {
        uint32 row = pos >> 5;
        uint64 val = this->extract2bits(pos);
        new_kmer.kmer_data[row] += val << (31 - (pos & 31)) * 2;
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
  virtual bool operator()(const VLMCKmer &a, const VLMCKmer &b) const = 0;
  [[nodiscard]] virtual VLMCKmer min_value() const {
    return VLMCKmer(0, 0, {});
  }
  [[nodiscard]] virtual VLMCKmer max_value() const {
    return VLMCKmer(MAX_K, 0, {});
  }
};

template <int MAX_K>
struct KMerComparator : public VirtualKMerComparator<MAX_K> {
  bool operator()(const VLMCKmer &a, const VLMCKmer &b) const override {
    return a < b;
  }
  [[nodiscard]] VLMCKmer min_value() const override {
    std::array<uint64, 4> vec{};
    return VLMCKmer(0, 0, vec);
  }
  [[nodiscard]] VLMCKmer max_value() const override {
    VLMCKmer max_kmer(MAX_K, (size_t)-1, {});
    for (int row = 0; row < max_kmer.n_rows; row++) {
      max_kmer.kmer_data[row] = 0xFFFFFFFFFFFFFFFF;
    }
    return max_kmer;
  }
};

template <int MAX_K>
struct ReverseKMerComparator : public VirtualKMerComparator<MAX_K> {
  bool operator()(const VLMCKmer &a, const VLMCKmer &b) const override {
    return a.reverse_less_than(b);
  }
  [[nodiscard]] VLMCKmer min_value() const override {
    VLMCKmer max_kmer(MAX_K, (size_t)-1, {});
    for (int row = 0; row < max_kmer.n_rows; row++) {
      max_kmer.kmer_data[row] = 0xFFFFFFFFFFFFFFFF;
    }
    return max_kmer;
  }
  [[nodiscard]] VLMCKmer max_value() const override {
    std::array<uint64, 4> vec{};
    return VLMCKmer(0, 0, vec);
  }
};

template <int MAX_K> struct KMerReverseKeyExtractor {
  typedef uint64 key_type;
  uint64 upper_bound = (1 - std::pow(4, MAX_K + 1)) / (1 - 4);

  [[nodiscard]] uint64 n_kmers_with_length(uint32_t len) const {
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

  [[nodiscard]] VLMCKmer min_value() const {
    VLMCKmer max_kmer(MAX_K, (size_t)-1, {});
    for (int row = 0; row < max_kmer.n_rows; row++) {
      max_kmer.kmer_data[row] = 0xFFFFFFFFFFFFFFFF;
    }
    return max_kmer;
  }
  [[nodiscard]] VLMCKmer max_value() const {
    std::array<uint64, 4> vec{};
    return VLMCKmer(0, 0, vec);
  }
};
