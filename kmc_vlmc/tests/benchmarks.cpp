#include <benchmark/benchmark.h>
#include <random>
#include <vector>

#include <kmc_file.h>

#include "../src/kmer.hpp"
#include "../src/support_pruning.hpp"

class KmerVLMCBenchmarks : public benchmark::Fixture {
public:
  void SetUp(const ::benchmark::State &state) override {
    auto status = kmer_database.OpenForListing("res");
  }
  CKMCFile kmer_database;

  VLMCTranslator kmer_api{12};
  VLMCKmer kmer{};
  uint64 counter;
};

BENCHMARK_F(KmerVLMCBenchmarks, FullKMCIteration)(benchmark::State &state) {
  size_t sum = 0.0;
  for (auto _ : state) {
    while (kmer_database.ReadNextKmer(kmer_api, counter)) {
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
      kmer = kmer_api.construct_vlmc_kmer();
      kmer.count = counter;
      sum += counter;
    }

    kmer_database.RestartListing();
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
      kmer = kmer_api.construct_vlmc_kmer();
      kmer.count = counter;
      benchmark::DoNotOptimize(
          VLMCKmer::get_first_differing_position(kmer_og, kmer));
      sum += counter;
    }

    kmer_database.RestartListing();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, TranslateKmers)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");

  size_t sum = 0.0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_translator.construct_vlmc_kmer());
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, DiffKmers)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  kmer_translator.from_string("ACTGACTGGTTT");
  VLMCKmer kmer_to = kmer_translator.construct_vlmc_kmer();

  size_t sum = 0.0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        VLMCKmer::get_first_differing_position(kmer_from, kmer_to));
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, CreatePrefixKmers)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  size_t sum = 0.0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(VLMCKmer::create_prefix_kmer(
        kmer_from, 4 + 1, 0, std::array<size_t, 4>{1, 2, 3, 4}));
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, WriteToSorter)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  OutOfCoreKmerContainer sorter{};

  size_t sum = 0.0;
  for (auto _ : state) {
    sorter.push(kmer_from);
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, OnlyKMCIteration)(benchmark::State &state) {
  size_t sum = 0.0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_database.ReadNextKmer(kmer_api, counter));
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, SupportPruningOutOfCore)
(benchmark::State &state) {
  OutOfCoreKmerContainer sorter{};

  size_t sum = 0.0;
  for (auto _ : state) {
    support_pruning<12>(kmer_database, sorter, 12,
                        [](int length, size_t count) -> bool { return true; });

    state.PauseTiming();
    kmer_database.RestartListing();
    state.ResumeTiming();
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, SupportPruningInCore)
(benchmark::State &state) {
  InCoreKmerContainer sorter{};

  size_t sum = 0.0;
  for (auto _ : state) {
    support_pruning<12>(kmer_database, sorter, 12,
                        [](int length, size_t count) -> bool { return true; });

    state.PauseTiming();
    kmer_database.RestartListing();
    state.ResumeTiming();
  }
}

BENCHMARK_MAIN();

//  increase_counts(counters, counter);
//
//  while (kmer_database.ReadNextKmer(kmer_api, counter)) {
//    kmer = kmer_api.construct_vlmc_kmer();
//    kmer.count = counter;
//    process_kmer(kmer, prev_kmer, counters, sorter, include_node, 1);
//
//    prev_kmer = kmer;
//  }
//  VLMCKmer epsilon(0, 0, {});
//  process_kmer(epsilon, prev_kmer, counters, sorter, include_node, 1);
//  output_root(counters[0], sorter);