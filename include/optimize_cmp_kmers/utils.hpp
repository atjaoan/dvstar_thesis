#pragma once

namespace utils {
  
void matrix_recursion(size_t start_index_left, size_t stop_index_left, size_t start_index_right, size_t stop_index_right, 
  const std::function<void(size_t &left, size_t &right)> &fun){
  auto diff_left = stop_index_left - start_index_left;
  auto diff_right = stop_index_right - start_index_right;
  if (diff_left == 1 && diff_right == 1){
    fun(start_index_left, start_index_right); 
  } else if (diff_right > diff_left){
    auto new_right_index = (stop_index_right + start_index_right) / 2;
    matrix_recursion(start_index_left, stop_index_left, start_index_right, new_right_index, fun);
    matrix_recursion(start_index_left, stop_index_left, new_right_index, stop_index_right, fun);
  } else {
    auto new_left_index = (stop_index_left + start_index_left) / 2;
    matrix_recursion(start_index_left, new_left_index, start_index_right, stop_index_right, fun);
    matrix_recursion(new_left_index, stop_index_left, start_index_right, stop_index_right, fun);
  }
}
}