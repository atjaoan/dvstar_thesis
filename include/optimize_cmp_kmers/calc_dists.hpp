#include <Eigen/Dense>
#include "cluster_container.hpp"
#include "vlmc_container.hpp"
#include "parallel.hpp"

namespace calculate {

using matrix_t  = Eigen::MatrixXd;
using vlmc_c    = container::VLMC_Container;
using cluster_c = container::Cluster_Container;   

void calculate_slice(size_t start_index, size_t stop_index, matrix_t &distances,
                     cluster_c &cluster_left, cluster_c &cluster_right,
                     const std::function<double(vlmc_c &, vlmc_c &)> &fun) {

  for (size_t i = start_index; i < stop_index; i++) {
    for (size_t j = 0; j < cluster_right.size(); j++) {
      distances(i, j) = fun(cluster_left.get(i), cluster_right.get(j));
    }
  }
}

matrix_t calculate_distances(
    cluster_c &cluster, std::function<double(vlmc_c &, vlmc_c &)> &distance_function,
    size_t requested_cores){

  matrix_t distances{cluster.size(), cluster.size()};

  auto fun = [&](size_t start_index, size_t stop_index) {
    calculate_slice(start_index, stop_index, distances,
                           cluster, cluster, distance_function);
  };
  //TODO use parallelize
  parallel::sequential(cluster.size(), fun, requested_cores);
  return distances; 
}

matrix_t calculate_distances(
    cluster_c &cluster_left, cluster_c &cluster_right,
    std::function<double(vlmc_c &, vlmc_c &)> &distance_function,
    size_t requested_cores){

      matrix_t distances{cluster_left.size(), cluster_right.size()};

      auto fun = [&](size_t start_index, size_t stop_index) {
      calculate_slice(start_index, stop_index, distances,
                           cluster_left, cluster_right, distance_function);
      };
      //TODO use parallelize
      parallel::sequential(cluster_left.size(), cluster_right.size(), fun, requested_cores);
      return distances; 
}
}