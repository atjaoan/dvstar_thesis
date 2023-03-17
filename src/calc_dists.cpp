#include "parser.hpp"
#include "get_cluster.hpp"
#include "calc_dists.hpp"

using matrix_t = Eigen::MatrixXd;

template <typename VC>
matrix_t calculate_cluster_distance(parser::cli_arguments arguments, const size_t nr_cores){
  auto distance_function = parser::parse_distance_function(arguments);
  auto cluster = cluster::get_cluster<VC>(arguments.first_VLMC_path, nr_cores, arguments.background_order, arguments.set_size);
  if (arguments.second_VLMC_path.empty()){
    return calculate::calculate_distances<VC>(cluster, distance_function, nr_cores);
  } else {
    auto cluster_to = cluster::get_cluster<VC>(arguments.second_VLMC_path, nr_cores, arguments.background_order, arguments.set_size);
    return calculate::calculate_distances<VC>(cluster, cluster_to, distance_function, nr_cores);
  }
}

matrix_t apply_container(parser::cli_arguments arguments, parser::VLMC_Rep vlmc_container, const size_t nr_cores){
  if (vlmc_container==parser::VLMC_Rep::vlmc_vector){
    return calculate_cluster_distance<container::VLMC_vector>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_indexing){
    return calculate_cluster_distance<container::VLMC_Indexing>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_sorted_vector){
    return calculate_cluster_distance<container::VLMC_sorted_vector>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_b_tree){
    return calculate_cluster_distance<container::VLMC_B_tree>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_hashmap){
    return calculate_cluster_distance<container::VLMC_hashmap>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_combo){
    return calculate_cluster_distance<container::VLMC_Combo>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_veb){
    return calculate_cluster_distance<container::VLMC_Veb>(arguments, nr_cores);
  }
}

matrix_t calculate_kmer_major(parser::cli_arguments arguments, const size_t nr_cores){
  //TODO use cores
  auto cluster = cluster::get_kmer_cluster(arguments.first_VLMC_path, arguments.background_order);
  if (arguments.second_VLMC_path.empty()){
    return calculate::calculate_distance_major(cluster, nr_cores);
  } else {
    auto cluster_to = cluster::get_kmer_cluster(arguments.second_VLMC_path);
    return calculate::calculate_distance_major(cluster, cluster_to, nr_cores);
  }
}

void print_matrix(matrix_t distance_matrix){
  for (size_t i = 0; i < distance_matrix.rows(); i++)
  {
    for (size_t j = 0; j < distance_matrix.cols(); j++)
    {
      std::cout << distance_matrix(i,j) << " ";
    }
    std::cout << std::endl;
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
  if(arguments.first_VLMC_path.empty()){
    std::cerr
        << "Error: A input path to .bintree files has to be given for comparison operation."
        << std::endl;
    return EXIT_FAILURE;
  }

  size_t nr_cores = parser::parse_dop(arguments.dop);

  if(arguments.mode==parser::Mode::compare){
    std::cout << "Running Normal implementation" << std::endl; 
    matrix_t distance_matrix = apply_container(arguments, arguments.vlmc, nr_cores);
    print_matrix(distance_matrix); 
  } else {
    std::cout << "Running Kmer-Major implemenation" << std::endl; 
    matrix_t distance_matrix = calculate_kmer_major(arguments, nr_cores);
    print_matrix(distance_matrix); 
  }

  return EXIT_SUCCESS;
}