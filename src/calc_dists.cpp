#include "parser.hpp"
#include "get_cluster.hpp"
#include "calc_dists.hpp"

using matrix_t = Eigen::MatrixXd;
using vlmc_c = container::Index_by_value;
using cluster_c = container::Cluster_Container;

void get_cluster_with_argument(container::Cluster_vector &cluster, std::filesystem::path path, parser::Vlmc_Rep vlmc_container){
  if (vlmc_container==parser::Vlmc_Rep::vlmc_vector){
    cluster::get_cluster<container::VLMC_vector>(path, cluster); 
  } else if (vlmc_container==parser::Vlmc_Rep::vlmc_multi_vector){
    cluster::get_cluster<container::Index_by_value>(path, cluster);
  } else if (vlmc_container==parser::Vlmc_Rep::vlmc_sorted_vector){
    cluster::get_cluster<container::VLMC_sorted_vector>(path, cluster);
  } else if (vlmc_container==parser::Vlmc_Rep::vlmc_b_tree){
    cluster::get_cluster<container::VLMC_B_tree>(path, cluster);
  } else if (vlmc_container==parser::Vlmc_Rep::vlmc_hashmap){
    cluster::get_cluster<container::VLMC_hashmap>(path, cluster);
  } else {
    std::cerr
      << "Error: The vlmc representation '" << vlmc_container << "' is not one that can be used :( "
      << std::endl;
    return EXIT_FAILURE;
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
  auto distance_function = parser::parse_distance_function(arguments.dist_fn, arguments.background_order);
  size_t nr_cores_to_use = parser::parse_dop(arguments.dop);
  
  if(arguments.second_VLMC_path.empty()){
    if (arguments.cluster == parser::Cluster_Rep::cluster_vector){
      container::Cluster_vector cluster{};
      get_cluster_with_argument(cluster, arguments.first_VLMC_path, arguments.vlmc); 
      matrix_t distance_matrix = calculate::calculate_distances(cluster, distance_function, nr_cores_to_use);
      
      for (size_t i = 0; i < distance_matrix.rows(); i++)
      {
        for (size_t j = 0; j < distance_matrix.cols(); j++)
        {
          std::cout << distance_matrix(i,j) << " ";
        }
        std::cout << std::endl;
      }
    } else {
      std::cerr
        << "Error: cluster representation is not one that can be used :( "
        << std::endl;
      return EXIT_FAILURE;
    }
  } else {
    if (arguments.cluster == parser::Cluster_Rep::cluster_vector){
      container::Cluster_vector left_cluster{};
      container::Cluster_vector right_cluster{};

      get_cluster_with_argument(left_cluster, arguments.first_VLMC_path, arguments.vlmc);
      get_cluster_with_argument(right_cluster, arguments.second_VLMC_path, arguments.vlmc);  

      matrix_t distance_matrix = calculate::calculate_distances(left_cluster, right_cluster, distance_function, nr_cores_to_use);

      for (size_t i = 0; i < distance_matrix.rows(); i++)
      {
        for (size_t j = 0; j < distance_matrix.cols(); j++)
        {
          std::cout << distance_matrix(i,j) << " ";
        }
        std::cout << std::endl;
      }
    } else {
      std::cerr
        << "Error: cluster representation is not one that can be used :( "
        << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}