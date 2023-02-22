#include "parser.hpp"
#include "get_cluster.hpp"
#include "calc_dists.hpp"

using matrix_t = Eigen::MatrixXd;

matrix_t calculate_cluster_distance(parser::cli_arguments arguments, parser::Vlmc_Rep vlmc_container, size_t nr_cores_to_use){
  if (vlmc_container==parser::Vlmc_Rep::vlmc_vector){
    auto distance_function = parser::parse_distance_function(arguments.dist_fn, arguments.background_order);
    auto cluster = cluster::get_cluster<container::VLMC_vector>(arguments.first_VLMC_path);
    if (arguments.second_VLMC_path.empty()){
      return calculate::calculate_distances<container::VLMC_vector>(cluster, distance_function, nr_cores_to_use);
    } else {
      auto cluster_to = cluster::get_cluster<container::VLMC_vector>(arguments.second_VLMC_path);
      return calculate::calculate_distances<container::VLMC_vector>(cluster, cluster_to, distance_function, nr_cores_to_use);
    }
  } else if (vlmc_container==parser::Vlmc_Rep::vlmc_multi_vector){
    auto distance_function = parser::parse_distance_function(arguments.dist_fn, arguments.background_order);
    auto cluster = cluster::get_cluster<container::Index_by_value>(arguments.first_VLMC_path);
    if (arguments.second_VLMC_path.empty()){
      return calculate::calculate_distances<container::Index_by_value>(cluster, distance_function, nr_cores_to_use);
    } else {
      auto cluster_to = cluster::get_cluster<container::Index_by_value>(arguments.second_VLMC_path);
      return calculate::calculate_distances<container::Index_by_value>(cluster, cluster_to, distance_function, nr_cores_to_use);
    }
  } else if (vlmc_container==parser::Vlmc_Rep::vlmc_sorted_vector){
    auto distance_function = parser::parse_distance_function(arguments.dist_fn, arguments.background_order);
    auto cluster = cluster::get_cluster<container::VLMC_sorted_vector>(arguments.first_VLMC_path);
    if (arguments.second_VLMC_path.empty()){
      return calculate::calculate_distances<container::VLMC_sorted_vector>(cluster, distance_function, nr_cores_to_use);
    } else {
      auto cluster_to = cluster::get_cluster<container::VLMC_sorted_vector>(arguments.second_VLMC_path);
      return calculate::calculate_distances<container::VLMC_sorted_vector>(cluster, cluster_to, distance_function, nr_cores_to_use);
    }
  } else if (vlmc_container==parser::Vlmc_Rep::vlmc_b_tree){
    auto distance_function = parser::parse_distance_function(arguments.dist_fn, arguments.background_order);
    auto cluster = cluster::get_cluster<container::VLMC_B_tree>(arguments.first_VLMC_path);
    if (arguments.second_VLMC_path.empty()){
      return calculate::calculate_distances<container::VLMC_B_tree>(cluster, distance_function, nr_cores_to_use);
    } else {
      auto cluster_to = cluster::get_cluster<container::VLMC_B_tree>(arguments.second_VLMC_path);
      return calculate::calculate_distances<container::VLMC_B_tree>(cluster, cluster_to, distance_function, nr_cores_to_use);
    }
  } else if (vlmc_container==parser::Vlmc_Rep::vlmc_hashmap){
    auto distance_function = parser::parse_distance_function(arguments.dist_fn, arguments.background_order);
    auto cluster = cluster::get_cluster<container::VLMC_hashmap>(arguments.first_VLMC_path);
    if (arguments.second_VLMC_path.empty()){
      return calculate::calculate_distances<container::VLMC_hashmap>(cluster, distance_function, nr_cores_to_use);
    } else {
      auto cluster_to = cluster::get_cluster<container::VLMC_hashmap>(arguments.second_VLMC_path);
      return calculate::calculate_distances<container::VLMC_hashmap>(cluster, cluster_to, distance_function, nr_cores_to_use);
    }
  }
}

int main(int argc, char *argv[]){
  CLI::App app{"Distance comparison of either one directory or between two different directories."};

  parser::cli_arguments arguments{};
  add_options(app, arguments);

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }
  if(arguments.mode != parser::Mode::compare){
    return EXIT_SUCCESS;
  }
  if(arguments.first_VLMC_path.empty()){
    std::cerr
        << "Error: A input path to .bintree files has to be given for comparison operation."
        << std::endl;
    return EXIT_FAILURE;
  }

  size_t nr_cores_to_use = parser::parse_dop(arguments.dop);

  matrix_t distance_matrix = calculate_cluster_distance(arguments, arguments.vlmc, nr_cores_to_use);
      
  for (size_t i = 0; i < distance_matrix.rows(); i++)
  {
    for (size_t j = 0; j < distance_matrix.cols(); j++)
    {
      std::cout << distance_matrix(i,j) << " ";
    }
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}