#include <benchmark/benchmark.h>
#include <random>
#include <vector>

#include <kmc_file.h>

#include "../src/build_vlmc.hpp"
#include "../src/kmer.hpp"
#include "../src/support_pruning.hpp"
#include "read_helper.hpp"

using namespace vlmc;

class KmerVLMCBenchmarks : public benchmark::Fixture {
  static void config_stxxl() {
    static bool have_configured = false;
    static std::string tmp_path = "test_tmp";
    if (!have_configured) {
      have_configured = true;
      std::filesystem::create_directories(tmp_path);
      configure_stxxl(tmp_path);
    }
  }

public:
  void SetUp(const ::benchmark::State &state) override {
    auto status = kmer_database.OpenForListing("res");
    config_stxxl();
  }

  void TearDown(const ::benchmark::State &state) override {
    kmer_database.Close();
  }

  CKMCFile kmer_database;

  VLMCTranslator kmer_api{12};
  VLMCKmer kmer{12, 0, {}};

  VLMCKmer prev_kmer = VLMCKmer{0, 0, {}};
  KMersPerLevel<uint64> counters{4 + 1, 12 + 1};
  const std::function<bool(int, size_t)> include_node =
      [](int length, size_t count) -> bool { return true; };

  const std::function<void(const VLMCKmer &, int, uint64,
                           std::vector<uint64> &)>
      output_func = [&](const VLMCKmer &prev_kmer_, int diff_pos, uint64 sum,
                        std::vector<uint64> &next_counts) {
        output_node(prev_kmer_, diff_pos, sum, next_counts, container);
      };

  vector_type container = std::make_shared<InCoreKmerContainer<>>();

  uint64 counter{};
};

BENCHMARK_F(KmerVLMCBenchmarks, TranslateKmers)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_translator.construct_vlmc_kmer());
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, AllocateKmer)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");

  for (auto _ : state) {
    benchmark::DoNotOptimize(VLMCKmer{12, 0, {}});
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, TranslateKmersWihoutAllocation)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");

  VLMCKmer kmer{12, 0, {}};
  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_translator.construct_vlmc_kmer(kmer));
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, DiffKmers)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  kmer_translator.from_string("ACTGACTGGTTT");
  VLMCKmer kmer_to = kmer_translator.construct_vlmc_kmer();

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        VLMCKmer::get_first_differing_position(kmer_from, kmer_to));
  }
}

static void PrefixKmers(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  auto diff_pos = state.range(0);

  for (auto _ : state) {
    benchmark::DoNotOptimize(VLMCKmer::create_prefix_kmer(
        kmer_from, diff_pos + 1, 0, std::array<uint64, 4>{1, 2, 3, 4}));
  }
}

BENCHMARK(PrefixKmers)->DenseRange(0, 11, 1);

static void PrefixKmersNoAllocation(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  VLMCKmer prefix_kmer{0, 0, {}};
  std::array<uint64, 4> counts{1, 2, 3, 4};

  auto diff_pos = state.range(0);

  for (auto _ : state) {
    benchmark::DoNotOptimize(VLMCKmer::create_prefix_kmer(
        kmer_from, diff_pos + 1, 0, counts, prefix_kmer));
  }
}

BENCHMARK(PrefixKmersNoAllocation)->DenseRange(0, 11, 1);

BENCHMARK_F(KmerVLMCBenchmarks, WriteToSorterOutOfCore)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  OutOfCoreKmerContainer<ReverseKMerComparator<31>> sorter{};

  for (auto _ : state) {
    sorter.push(kmer_from);
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, WriteToSorterInCore)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  InCoreKmerContainer<ReverseKMerComparator<31>> sorter{};

  for (auto _ : state) {
    sorter.push(kmer_from);
  }
}


BENCHMARK_F(KmerVLMCBenchmarks, OnlyKMCIteration)(benchmark::State &state) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_database.ReadNextKmer(kmer_api, counter));
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, PrefixIndex)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer = kmer_translator.construct_vlmc_kmer();

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer.get_prefix_index(3));
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, CompareReverse)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ATCGATCGACTT");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  kmer_translator.from_string("CGACGATCAGCA");
  VLMCKmer kmer_to = kmer_translator.construct_vlmc_kmer();

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_from.reverse_less_than(kmer_to));
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, Compare)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ATCGATCGACTT");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  kmer_translator.from_string("CGACGATCAGCA");
  VLMCKmer kmer_to = kmer_translator.construct_vlmc_kmer();

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_from < kmer_to);
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, CharPos)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ATCGATCGACTT");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_from.char_pos(5));
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, Extract2Bits)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ATCGATCGACTT");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_from.extract2bits(5));
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, IncreaseCounters)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ATCGATCGACTT");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();
  kmer_from.count = 50;
  KMersPerLevel<uint64> counters{4 + 1, 12 + 1};

  for (auto _ : state) {
    increase_counts(counters, kmer_from);
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, SupportPruningSplitTasksInCore)
(benchmark::State &state) {
  std::shared_ptr<KmerContainer<ReverseKMerComparator<31>>> sorter =
      std::make_shared<InCoreKmerContainer<ReverseKMerComparator<31>>>();

  int prefix_length = 3;
  int n_tasks = VLMCKmer::n_kmers_with_length(prefix_length);

  std::vector<std::shared_ptr<KmerContainer<KMerComparator<31>>>> tasks(
      n_tasks);
  for (int i = 0; i < n_tasks; i++) {
    tasks[i] = std::make_shared<InCoreKmerContainer<KMerComparator<31>>>();
  }

  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      kmer = kmer_api.construct_vlmc_kmer();
      kmer.count = counter;

      size_t prefix_idx = kmer.get_prefix_index(prefix_length);

      tasks[prefix_idx]->push(kmer);
    }

    state.PauseTiming();
    kmer_database.RestartListing();
    state.ResumeTiming();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, SupportPruningSplitTasksOutOfCore)
(benchmark::State &state) {
  int prefix_length = 3;
  int n_tasks = VLMCKmer::n_kmers_with_length(prefix_length);

  std::vector<std::shared_ptr<KmerContainer<KMerComparator<31>>>> tasks(
      n_tasks);
  for (int i = 0; i < n_tasks; i++) {
    tasks[i] = std::make_shared<OutOfCoreKmerContainer<KMerComparator<31>>>();
  }

  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      kmer = kmer_api.construct_vlmc_kmer();
      kmer.count = counter;

      size_t prefix_idx = kmer.get_prefix_index(prefix_length);

      tasks[prefix_idx]->push(kmer);
    }

    state.PauseTiming();
    kmer_database.RestartListing();
    state.ResumeTiming();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, SupportPruningOutOfCore)
(benchmark::State &state) {
  std::shared_ptr<KmerContainer<ReverseKMerComparator<31>>> sorter =
      std::make_shared<OutOfCoreKmerContainer<ReverseKMerComparator<31>>>();

  for (auto _ : state) {
    support_pruning<12>(
        kmer_database, sorter, 12,
        [](int length, size_t count) -> bool { return true; }, "external");

    state.PauseTiming();
    kmer_database.RestartListing();
    state.ResumeTiming();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, SupportPruningInCore)
(benchmark::State &state) {
  std::shared_ptr<KmerContainer<ReverseKMerComparator<31>>> sorter =
      std::make_shared<InCoreKmerContainer<ReverseKMerComparator<31>>>();

  for (auto _ : state) {
    support_pruning<12>(
        kmer_database, sorter, 12,
        [](int length, size_t count) -> bool { return true; }, "internal");

    state.PauseTiming();
    kmer_database.RestartListing();
    state.ResumeTiming();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, SequentialSupportPruningOutOfCore)
(benchmark::State &state) {
  std::shared_ptr<KmerContainer<ReverseKMerComparator<31>>> sorter =
      std::make_shared<OutOfCoreKmerContainer<ReverseKMerComparator<31>>>();

  for (auto _ : state) {
    sequential_support_pruning<12>(
        kmer_database, sorter, 12,
        [](int length, size_t count) -> bool { return true; });

    state.PauseTiming();
    kmer_database.RestartListing();
    state.ResumeTiming();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, SequentialSupportPruningInCore)
(benchmark::State &state) {
  std::shared_ptr<KmerContainer<ReverseKMerComparator<31>>> sorter =
      std::make_shared<InCoreKmerContainer<ReverseKMerComparator<31>>>();

  for (auto _ : state) {
    sequential_support_pruning<12>(
        kmer_database, sorter, 12,
        [](int length, size_t count) -> bool { return true; });

    state.PauseTiming();
    kmer_database.RestartListing();
    state.ResumeTiming();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, FullSupportPruning)
(benchmark::State &state) {
  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      kmer_api.construct_vlmc_kmer(kmer);
      kmer.count = counter;
      if (prev_kmer.length != 0) {
        process_kmer(kmer, prev_kmer, counters, output_func, include_node, 1);
      }
      increase_counts(counters, kmer);
      prev_kmer = kmer;
    }

    //    kmer_database.RestartListing();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, ManualSupportPruningWithoutOutputWithPrefix)
(benchmark::State &state) {
  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      kmer_api.construct_vlmc_kmer(kmer);
      kmer.count = counter;

      if (prev_kmer.length != 0) {
        VLMCKmer output_kmer{0, 0, {}};

        int32 diff_pos = 0;
        if (kmer.length == 0) {
          // Output the longest prefix of all the kmers in the current set
          diff_pos = 0;
        } else {
          diff_pos = VLMCKmer::get_first_differing_position(kmer, prev_kmer);
        }

        if (diff_pos != -1) {

          // Save the counts of the k-mers that are to be outputted:
          for (int32 i = prev_kmer.length - 1; i >= diff_pos; i--) {
            // For pseudo-counts, add 4 to counter
            // And the length has to be at least one shorter than the prev_kmer
            // to allow variable k iterations
            auto &next_counts = counters[diff_pos + 1];

            uint64 sum = next_counts[1] + next_counts[2] + next_counts[3] +
                         next_counts[4];
            auto prefix_kmer = VLMCKmer::create_prefix_kmer(
                prev_kmer, diff_pos + 1, sum,
                std::array<uint64, 4>{next_counts[1], next_counts[2],
                                      next_counts[3], next_counts[4]});

            auto char_idx = prev_kmer.char_pos(i);
            uchar offset_char_idx = char_idx + 1;
            counters[i][offset_char_idx] = sum;

            // Reset counters for potentially outputted kmer
            next_counts[1] = next_counts[2] = next_counts[3] = next_counts[4] =
                0;
          }
        }
      }

      increase_counts(counters, kmer);
      prev_kmer = kmer;
    }

    //    kmer_database.RestartListing();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, ManualSupportPruningWithoutOutput)
(benchmark::State &state) {
  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      kmer_api.construct_vlmc_kmer(kmer);
      kmer.count = counter;

      if (prev_kmer.length != 0) {
        VLMCKmer output_kmer{0, 0, {}};

        int32 diff_pos = 0;
        if (kmer.length == 0) {
          // Output the longest prefix of all the kmers in the current set
          diff_pos = 0;
        } else {
          diff_pos = VLMCKmer::get_first_differing_position(kmer, prev_kmer);
        }

        if (diff_pos != -1) {

          // Save the counts of the k-mers that are to be outputted:
          for (int32 i = prev_kmer.length - 1; i >= diff_pos; i--) {
            // For pseudo-counts, add 4 to counter
            // And the length has to be at least one shorter than the prev_kmer
            // to allow variable k iterations
            auto &next_counts = counters[diff_pos + 1];

            uint64 sum = next_counts[1] + next_counts[2] + next_counts[3] +
                         next_counts[4];

            auto char_idx = prev_kmer.char_pos(i);
            uchar offset_char_idx = char_idx + 1;
            counters[i][offset_char_idx] = sum;

            // Reset counters for potentially outputted kmer
            next_counts[1] = next_counts[2] = next_counts[3] = next_counts[4] =
                0;
          }
        }
      }

      increase_counts(counters, kmer);
      prev_kmer = kmer;
    }

    //    kmer_database.RestartListing();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, ManualSupportPruningWithoutOutputAndReset)
(benchmark::State &state) {
  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      kmer_api.construct_vlmc_kmer(kmer);
      kmer.count = counter;

      if (prev_kmer.length != 0) {
        VLMCKmer output_kmer{0, 0, {}};

        int32 diff_pos = 0;
        if (kmer.length == 0) {
          // Output the longest prefix of all the kmers in the current set
          diff_pos = 0;
        } else {
          diff_pos = VLMCKmer::get_first_differing_position(kmer, prev_kmer);
        }

        if (diff_pos != -1) {

          // Save the counts of the k-mers that are to be outputted:
          // Save the counts of the k-mers that are to be outputted:
          for (int32 i = prev_kmer.length - 1; i >= diff_pos; i--) {
            // For pseudo-counts, add 4 to counter
            // And the length has to be at least one shorter than the prev_kmer
            // to allow variable k iterations
            uint64 sum = counters[i + 1][1] + counters[i + 1][2] +
                         counters[i + 1][3] + counters[i + 1][4];

            auto char_idx = prev_kmer.char_pos(i);
            uchar offset_char_idx = char_idx + 1;
            counters[i][offset_char_idx] = sum;
          }
        }
      }

      increase_counts(counters, kmer);
      prev_kmer = kmer;
    }

    //    kmer_database.RestartListing();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks,
            ManualSupportPruningWithoutOutputAndResetAndLoop)
(benchmark::State &state) {
  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      kmer_api.construct_vlmc_kmer(kmer);
      kmer.count = counter;

      if (prev_kmer.length != 0) {
        VLMCKmer output_kmer{0, 0, {}};

        int32 diff_pos = 0;
        if (kmer.length == 0) {
          // Output the longest prefix of all the kmers in the current set
          diff_pos = 0;
        } else {
          diff_pos = VLMCKmer::get_first_differing_position(kmer, prev_kmer);
        }

        if (diff_pos != -1) {

          for (int32 i = prev_kmer.length - 1; i >= diff_pos; i--) {
            // For pseudo-counts, add 4 to counter
            // And the length has to be at least one shorter than the prev_kmer
            // to allow variable k iterations
            uint64 sum;
            benchmark::DoNotOptimize(sum = 0);
          }
        }
      }

      increase_counts(counters, kmer);
      prev_kmer = kmer;
    }

    //    kmer_database.RestartListing();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks,
            ManualSupportPruningWithoutOutputAndResetAndLoopAndSaveAndDiffPos)
(benchmark::State &state) {
  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      kmer_api.construct_vlmc_kmer(kmer);
      kmer.count = counter;

      if (prev_kmer.length != 0) {
        VLMCKmer output_kmer{0, 0, {}};

        int32 diff_pos = 0;
        if (kmer.length == 0) {
          // Output the longest prefix of all the kmers in the current set
          benchmark::DoNotOptimize(diff_pos = 0);
        } else {
        }
      }

      increase_counts(counters, kmer);
      prev_kmer = kmer;
    }

    //    kmer_database.RestartListing();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, KMCIterationWithDiff)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer_og = kmer_translator.construct_vlmc_kmer();

  size_t sum = 0.0;
  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      kmer_api.construct_vlmc_kmer(kmer);
      kmer.count = counter;

      benchmark::DoNotOptimize(
          VLMCKmer::get_first_differing_position(kmer_og, kmer));
      sum += counter;
    }

    kmer_database.RestartListing();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, KMCIterationWithConversion)
(benchmark::State &state) {
  size_t sum = 0.0;
  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      kmer_api.construct_vlmc_kmer(kmer);
      kmer.count = counter;
      sum += counter;
    }

    kmer_database.RestartListing();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, FullKMCIteration)(benchmark::State &state) {
  size_t sum = 0.0;
  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
      sum += counter;
    }

    kmer_database.RestartListing();
  }
}

BENCHMARK_DEFINE_F(KmerVLMCBenchmarks, SupportPruningInnerLoopDiffNoOutput)
(benchmark::State &state) {
  VLMCKmer prefix_kmer{0, 0, {}};

  std::string kmer_string{"AAAAAAAAAAAA"};
  kmer = create_kmer(kmer_string);

  auto diff_index = state.range(0) - 1;
  kmer_string[diff_index] = 'C';
  prev_kmer = create_kmer(kmer_string);
  container->clear();

  for (auto _ : state) {
    int32 diff_pos = VLMCKmer::get_first_differing_position(kmer, prev_kmer);
    // Output kmers, longest to shortest to get counts right
    for (int i = prev_kmer.length - 1; i >= diff_pos; i--) {
      // For pseudo-counts, add 4 to counter
      // And the length has to be at least one shorter than the prev_kmer
      // to allow variable k iterations
      auto &next_counts = counters[i + 1];

      uint64 sum =
          next_counts[1] + next_counts[2] + next_counts[3] + next_counts[4];

      if (i < prev_kmer.length - 1 && include_node(i + 1, sum + 4)) {
        VLMCKmer::create_prefix_kmer(prev_kmer, i + 1, sum,
                                     std::array<uint64, 4>{1, 2, 3, 4},
                                     prefix_kmer);
        container->push(prefix_kmer);
      }

      auto char_idx = prev_kmer.char_pos(i);
      counters[i][char_idx + 1] = sum;

      // Reset counters for potentially outputted kmer
      counters[i][0] = 0;
      next_counts[1] = next_counts[2] = next_counts[3] = next_counts[4] = 0;
    }
  }
}

BENCHMARK_REGISTER_F(KmerVLMCBenchmarks, SupportPruningInnerLoopDiffNoOutput)
    ->DenseRange(1, 12, 1);

BENCHMARK_DEFINE_F(KmerVLMCBenchmarks, SupportPruningInnerLoopDiff)
(benchmark::State &state) {
  std::string kmer_string{"AAAAAAAAAAAA"};
  kmer = create_kmer(kmer_string);

  auto diff_index = state.range(0) - 1;
  kmer_string[diff_index] = 'C';
  prev_kmer = create_kmer(kmer_string);
  container->clear();

  for (auto _ : state) {

    int32 diff_pos = VLMCKmer::get_first_differing_position(kmer, prev_kmer);
    // Output kmers, longest to shortest to get counts right
    for (int i = prev_kmer.length - 1; i >= diff_pos; i--) {
      // For pseudo-counts, add 4 to counter
      // And the length has to be at least one shorter than the prev_kmer
      // to allow variable k iterations
      auto &next_counts = counters[i + 1];

      uint64 sum =
          next_counts[1] + next_counts[2] + next_counts[3] + next_counts[4];

      if (i < prev_kmer.length - 1 && include_node(i + 1, sum + 4)) {
        output_node(prev_kmer, i, sum, next_counts, container);
      }

      auto char_idx = prev_kmer.char_pos(i);
      counters[i][char_idx + 1] = sum;

      // Reset counters for potentially outputted kmer
      counters[i][0] = 0;
      next_counts[1] = next_counts[2] = next_counts[3] = next_counts[4] = 0;
    }
  }
}

BENCHMARK_REGISTER_F(KmerVLMCBenchmarks, SupportPruningInnerLoopDiff)
    ->DenseRange(1, 12, 1);

BENCHMARK_MAIN();
