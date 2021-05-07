#include <benchmark/benchmark.h>
#include <random>
#include <vector>

#include <kmc_file.h>

#include "../src/build_vlmc.hpp"
#include "../src/kmer.hpp"
#include "../src/support_pruning.hpp"

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
  VLMCKmer kmer{};
  uint64 counter{};
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
        kmer_from, 4 + 1, 0, std::array<uint64, 4>{1, 2, 3, 4}));
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, WriteToSorter)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  OutOfCoreKmerContainer<ReverseKMerComparator<31>> sorter{};

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

BENCHMARK_F(KmerVLMCBenchmarks, PrefixIndex)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ACTGACTGACTG");
  VLMCKmer kmer = kmer_translator.construct_vlmc_kmer();

  size_t sum = 0.0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer.get_prefix_index(3));
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

  size_t sum = 0.0;
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

  size_t sum = 0.0;
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

BENCHMARK_F(KmerVLMCBenchmarks, CompareReverse)
(benchmark::State &state) {
  VLMCTranslator kmer_translator{12};
  kmer_translator.from_string("ATCGATCGACTT");
  VLMCKmer kmer_from = kmer_translator.construct_vlmc_kmer();

  kmer_translator.from_string("CGACGATCAGCA");
  VLMCKmer kmer_to = kmer_translator.construct_vlmc_kmer();

  size_t sum = 0.0;
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

  size_t sum = 0.0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_from < kmer_to);
  }
}

BENCHMARK_F(KmerVLMCBenchmarks, SupportPruningOutOfCore)
(benchmark::State &state) {
  std::shared_ptr<KmerContainer<ReverseKMerComparator<31>>> sorter =
      std::make_shared<OutOfCoreKmerContainer<ReverseKMerComparator<31>>>();

  size_t sum = 0.0;
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

  size_t sum = 0.0;
  for (auto _ : state) {
    support_pruning<12>(
        kmer_database, sorter, 12,
        [](int length, size_t count) -> bool { return true; }, "internal");

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