#pragma once

#include <stdlib.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

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

void recursive_get_triangle_coords(std::vector<std::array<int, 6>> &triangle_coords, 
  int x1, int y1, int x2, int y2, int x3, int y3, const int remaining_cores){
  if (remaining_cores == 1) {
    triangle_coords.push_back({x1, y1, x2, y2, x3, y3}); 
    return; 
  }
  auto one_to_three  = x1 == x3 || y1 == y3;
  auto two_to_three  = x2 == x3 || y2 == y3;
  auto one_to_two    = x1 == x2 || y1 == y2; 

  if (two_to_three && one_to_two){
    auto new_x = (x3 + x1 + 1) / 2;
    auto new_y = (y3 + y1 + 1) / 2;
    if (remaining_cores <= 2) {
      triangle_coords.push_back({x1, y1, x2, y2, new_x, new_y});
      triangle_coords.push_back({new_x, new_y, x2, y2, x3, y3});
    } else {
      auto cores = remaining_cores - (remaining_cores / 2); 
      recursive_get_triangle_coords(triangle_coords, x1, y1, x2, y2, new_x, new_y, remaining_cores / 2);
      recursive_get_triangle_coords(triangle_coords, new_x, new_y, x2, y2, x3, y3, cores);
    }
  } else if ((two_to_three && not one_to_three) || (one_to_three && one_to_two)) {
    auto new_x = (x3 + x2 + 1) / 2;
    auto new_y = (y3 + y2 + 1) / 2; 
    if (remaining_cores <= 2) {
      triangle_coords.push_back({x1, y1, x2, y2, new_x, new_y});
      triangle_coords.push_back({x1, y1, new_x, new_y, x3, y3});
    } else {
      auto cores = remaining_cores - (remaining_cores / 2);
      recursive_get_triangle_coords(triangle_coords, x1, y1, x2, y2, new_x, new_y, remaining_cores / 2);
      recursive_get_triangle_coords(triangle_coords, x1, y1, new_x, new_y, x3, y3, cores);
    }
  } else {
    auto new_x = (x2 + x1 + 1) / 2; 
    auto new_y = (y2 + y1 + 1) / 2; 
    if (remaining_cores <= 2) {
      triangle_coords.push_back({x1, y1, new_x, new_y, x3, y3});
      triangle_coords.push_back({new_x, new_y, x2, y2, x3, y3});
    } else {
      auto cores = remaining_cores - (remaining_cores / 2);  
      recursive_get_triangle_coords(triangle_coords, x1, y1, new_x, new_y, x3, y3, remaining_cores / 2);
      recursive_get_triangle_coords(triangle_coords, new_x, new_y, x2, y2, x3, y3, cores);
    }
  }
}

void parallelize_triangle(size_t size, const std::function<void(int, int, int, int, int, int)> &fun, const size_t requested_cores) {
  std::vector<std::thread> threads{};
  std::vector<std::array<int, 6>> triangle_coords; 
  recursive_get_triangle_coords(triangle_coords, 0, 0, 0, size, size, size, requested_cores);

  for (auto &coords : triangle_coords) {
    threads.emplace_back(fun, coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
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

void pool_parallel(size_t size_left, const std::function<void(size_t)> &fun, BS::thread_pool& pool) {
  for (auto left = 0; left < size_left; left++) {
    std::future<void> fut = pool.submit(fun, left);
  }

  pool.wait_for_tasks();
}

void parallelize_kmer_major(size_t size, const std::function<void(size_t, size_t, size_t)> &fun, BS::thread_pool& pool) {
  std::vector<std::tuple<size_t, size_t>> bounds{};
  float values_per_thread = std::floor(std::sqrt(size));

  auto start_index = 0; 
  while (start_index < size) {
    if (start_index + values_per_thread > size) {
      bounds.emplace_back(start_index, size); 
      break; 
    }
    bounds.emplace_back(start_index, start_index + values_per_thread);
    start_index += values_per_thread; 
  }

  int idx = 0; 
  for (auto &[start_index, stop_index] : bounds) {
    std::future<void> fut = pool.submit(fun, start_index, stop_index, idx);
    idx++; 
  }
  pool.wait_for_tasks();
}
}