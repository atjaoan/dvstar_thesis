#include <Eigen/Dense>
#include "cluster_container.hpp"
#include "vlmc_container.hpp"
#include "parallel.hpp"

namespace calculate {

using matrix_t  = Eigen::MatrixXd;
using vlmc_c    = container::VLMC_Container;
using cluster_c = container::Cluster_Container;   

void calculate_slice(size_t start_index, size_t stop_index, matrix_t &distances,
                     cluster_c &trees, cluster_c &trees_to,
                     const std::function<float(vlmc_c &, vlmc_c &)> &fun) {

  for (size_t i = start_index; i < stop_index; i++) {
    for (size_t j = 0; j < trees_to.size(); j++) {
      distances(i, j) = fun(trees.get(i), trees_to.get(i));
    }
  }
}

matrix_t calculate_distances(
    cluster_c &trees, std::function<float(vlmc_c &, vlmc_c &)> &distance_function,
    size_t requested_cores){

  matrix_t distances{trees.size(), trees.size()};

  auto fun = [&](size_t start_index, size_t stop_index) {
    calculate_slice(start_index, stop_index, distances,
                           trees, trees, distance_function);
  };
  //TODO use parallelize
  parallel::sequential(trees.size(), fun, requested_cores);
  return distances; 
}

matrix_t calculate_distances(
    cluster_c &left_trees, cluster_c &right_trees,
    std::function<float(vlmc_c &, vlmc_c &)> &distance_function,
    size_t requested_cores){

      matrix_t distances{left_trees.size(), right_trees.size()};

      auto fun = [&](size_t start_index, size_t stop_index) {
      calculate_slice(start_index, stop_index, distances,
                           left_trees, right_trees, distance_function);
      };
      //TODO use parallelize
      parallel::sequential(left_trees.size(), right_trees.size(), fun, requested_cores);
      return distances; 
}
}