#pragma once

#include <stdlib.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

namespace parallel {

std::vector<std::tuple<size_t, size_t>> get_bounds(size_t size, const size_t requested_cores) {
  const size_t processor_count = std::thread::hardware_concurrency();
  std::vector<std::tuple<size_t, size_t>> bounds_per_thread{};
  float values_per_thread = float(size) / float(processor_count);

  auto limit = std::min(processor_count, size);
  for (size_t i = 0; i < limit; i++) {
    size_t start_index = std::floor(values_per_thread * i);
    size_t stop_index = std::floor(values_per_thread * (i + 1.0));
    if (i == (limit - 1)) {
      stop_index = size;
    }
    bounds_per_thread.emplace_back(start_index, stop_index);
  }

  return bounds_per_thread;
}

void parallelize(size_t size, const std::function<void(size_t, size_t)> &fun, const size_t requested_cores) {
  std::vector<std::thread> threads{};

  auto bounds = get_bounds(size, requested_cores);
  for (auto &[start_index, stop_index] : bounds) {
    threads.emplace_back(fun, start_index, stop_index);
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void parallelize(size_t size_left, size_t size_right, const std::function<void(size_t, size_t)> &fun, const size_t requested_cores) {
  std::vector<std::thread> threads{};
  // TODO only uses one size...
  auto bounds = get_bounds(size_left, requested_cores);
  for (auto &[start_index, stop_index] : bounds) {
    threads.emplace_back(fun, start_index, stop_index);
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void sequential(size_t size, const std::function<void(size_t, size_t)> &fun, const size_t requested_cores) {
  std::thread t {fun, size, size};
  t.join();
}

void sequential(size_t size_left, size_t size_right, const std::function<void(size_t, size_t)> &fun, const size_t requested_cores) {
  std::thread t {fun, size_left, size_right};
  t.join();
}

}