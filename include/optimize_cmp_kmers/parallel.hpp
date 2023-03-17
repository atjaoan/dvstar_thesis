#pragma once

#include <stdlib.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>
#include "BS_thread_pool.hpp"

#include "BS_thread_pool.hpp"

namespace parallel {

std::vector<std::tuple<size_t, size_t>> get_x_bounds(size_t size, const size_t requested_cores) {
  const size_t processor_count = std::thread::hardware_concurrency();
  size_t used_cores = 1;
  if(requested_cores > size){
    used_cores = size;
  } else if(requested_cores <= processor_count){
      used_cores = requested_cores;
  } else {
    used_cores = processor_count;
  }
  std::vector<std::tuple<size_t, size_t>> bounds_per_thread{};
  float values_per_thread = float(size) / float(used_cores);

  auto limit = std::min(used_cores, size);
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

std::vector<std::tuple<size_t, size_t>> get_exp_x_bounds(size_t size, const size_t requested_cores) {
  const size_t processor_count = std::thread::hardware_concurrency();
  size_t used_cores = 1;

  if(requested_cores <= processor_count) used_cores = requested_cores;
  else used_cores = processor_count;

  std::vector<std::tuple<size_t, size_t>> bounds_per_thread{};

  int i = size;
  for (; i > 20 && used_cores > 1; i = (int)std::floor(i*0.2)) {
    size_t start_index = (int)std::floor(i*0.2);
    size_t stop_index = i;
    bounds_per_thread.emplace_back(start_index, stop_index);
    used_cores--;
  }
  bounds_per_thread.emplace_back(0, i);
  return bounds_per_thread;
}

void parallelize(size_t size, const std::function<void(size_t, size_t)> &fun, const size_t requested_cores) {
  std::vector<std::thread> threads{};

  auto bounds = get_x_bounds(size, requested_cores);
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
  auto bounds = get_x_bounds(size_left, requested_cores);
  for (auto &[start_index, stop_index] : bounds) {
    threads.emplace_back(fun, start_index, stop_index);
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void parallelize_exp_halving(size_t size_left, const std::function<void(size_t, size_t)> &fun, const size_t requested_cores) {
  std::vector<std::thread> threads{};
  auto bounds = get_exp_x_bounds(size_left, requested_cores);
  for (auto &[start_index, stop_index] : bounds) {
    threads.emplace_back(fun, start_index, stop_index);
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void parallelize_exp_halving(size_t size_left, size_t size_right, const std::function<void(size_t, size_t)> &fun, const size_t requested_cores) {
  std::vector<std::thread> threads{};
  auto bounds = get_exp_x_bounds(size_left, requested_cores);
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
  std::thread t {fun, 0, size};
  t.join();
}

void sequential(size_t size_left, size_t size_right, const std::function<void(size_t, size_t)> &fun, const size_t requested_cores) {
  std::thread t {fun, 0, size_left};
  t.join();
}

void pool_parallelize(size_t size, const std::function<void(size_t, size_t)> &fun, const size_t requested_cores, BS::thread_pool& pool) {
  

  auto bounds = get_x_bounds(size, requested_cores);
  for (auto &[start_index, stop_index] : bounds) {
    std::future<void> fut = pool.submit(fun, start_index, stop_index);
  }
  pool.wait_for_tasks();
}

}