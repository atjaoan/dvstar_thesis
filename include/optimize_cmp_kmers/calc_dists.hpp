#include <Eigen/Dense>
#include "cluster_container.hpp"
#include "vlmc_container.hpp"
#include "parallel.hpp"

namespace calculate{

using matrix_t  = Eigen::MatrixXd;
using vlmc_c    = container::VLMC_Container;
using cluster_c = container::Cluster_Container;   

matrix_t calculate_distances(cluster_c &trees, std::function<float(vlmc_c &, vlmc_c &)> &distance_function){
  matrix_t distances{trees.size(), trees.size()};

  auto fun = [&](size_t start_index, size_t stop_index) {
    calculate_vector_slice(start_index, stop_index, std::ref(distances),
                           std::ref(trees), distance_function);
  };
  parallel::parallelize(trees.size(), fun);
  return distances; 
}
    
}