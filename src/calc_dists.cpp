#pragma once 

#include <highfive/H5File.hpp>

#include "parser.hpp"
#include "get_cluster.hpp"
#include "calc_dists.hpp"
#include "global_aliases.hpp"
#include "utils.hpp"

matrix_t calculate_kmer_major(parser::cli_arguments arguments, const size_t nr_cores) {
  size_t use_cores = nr_cores;
  size_t max_cores = std::thread::hardware_concurrency();
  if (max_cores < nr_cores) {
    use_cores = max_cores;
  }

  auto cluster = get_cluster::get_kmer_cluster(arguments.first_VLMC_path, use_cores, arguments.background_order, arguments.set_size);
  if (arguments.second_VLMC_path.empty()) {
    std::cout << "Calculating distances for single cluster." << std::endl;
    return calc_dist::calculate_distance_major(cluster, cluster, use_cores);
  }
  auto cluster_to = get_cluster::get_kmer_cluster(arguments.second_VLMC_path, use_cores, arguments.background_order, arguments.set_size);
  std::cout << "Calculating distances." << std::endl;
  return calc_dist::calculate_distance_major(cluster, cluster_to, use_cores);
}

template <typename VC>
matrix_t calculate_cluster_distance(parser::cli_arguments arguments, const size_t nr_cores) {
  auto cluster = get_cluster::get_cluster<VC>(arguments.first_VLMC_path, nr_cores, arguments.background_order, arguments.set_size);
  if (arguments.second_VLMC_path.empty()) {
    std::cout << "Calculating distances for single cluster of size " << cluster.size() << std::endl;
    return calc_dist::calculate_distances<VC>(cluster, nr_cores);
  }
  auto cluster_to = get_cluster::get_cluster<VC>(arguments.second_VLMC_path, nr_cores, arguments.background_order, arguments.set_size);
  std::cout << "Calculating distances matrix of size " << cluster.size() << "x" << cluster_to.size() << std::endl;
  return calc_dist::calculate_distances<VC>(cluster, cluster_to, nr_cores);
}

matrix_t apply_container(parser::cli_arguments arguments, parser::VLMC_Rep vlmc_container, const size_t nr_cores) {
  if (vlmc_container == parser::VLMC_Rep::vlmc_sorted_vector) {
    return calculate_cluster_distance<vlmc_container::VLMC_sorted_vector>(arguments, nr_cores);
  }
  else if (vlmc_container == parser::VLMC_Rep::vlmc_b_tree) {
    return calculate_cluster_distance<vlmc_container::VLMC_B_tree>(arguments, nr_cores);
  }
  else if (vlmc_container == parser::VLMC_Rep::vlmc_hashmap) {
    return calculate_cluster_distance<vlmc_container::VLMC_hashmap>(arguments, nr_cores);
  }
  else if (vlmc_container == parser::VLMC_Rep::vlmc_veb) {
    return calculate_cluster_distance<vlmc_container::VLMC_Veb>(arguments, nr_cores);
  }
  else if (vlmc_container == parser::VLMC_Rep::vlmc_ey) {
    return calculate_cluster_distance<vlmc_container::VLMC_Eytzinger>(arguments, nr_cores);
  }
  else if (vlmc_container == parser::VLMC_Rep::vlmc_sorted_search) {
    return calculate_cluster_distance<vlmc_container::VLMC_sorted_search>(arguments, nr_cores);
  }
  else if (vlmc_container == parser::VLMC_Rep::vlmc_kmer_major) {
    return calculate_kmer_major(arguments, nr_cores);
  }
}

int main(int argc, char* argv[]) {
  CLI::App app{"Distance comparison of either one or between two directories of VLMCs."};

  parser::cli_arguments arguments{};
  add_options(app, arguments);

  try {
    app.parse(argc, argv);
  }
  catch (const CLI::ParseError& e) {
    return app.exit(e);
  }
  if (arguments.first_VLMC_path.empty()) {
    std::cerr
      << "Error: A input path to .bintree files has to be given for comparison operation."
      << std::endl;
    return EXIT_FAILURE;
  }

  size_t nr_cores = parser::parse_dop(arguments.dop);

  matrix_t distance_matrix = apply_container(arguments, arguments.vlmc, nr_cores);

  if (arguments.out_path.empty()) {
    // utils::print_matrix(distance_matrix);
  }
  else if (arguments.out_path.extension() == ".h5" ||
    arguments.out_path.extension() == ".hdf5") {
    HighFive::File file{arguments.out_path, HighFive::File::OpenOrCreate};

    if (!file.exist("distances")) {
      file.createGroup("distances");
    }
    auto distance_group = file.getGroup("distances");

    if (!distance_group.exist("distances")) {
      std::vector<size_t> dims{distance_matrix.rows(), distance_matrix.cols()};
      distance_group.createDataSet<double>("distances",
        HighFive::DataSpace(dims));
    }

    auto distance_data_set = distance_group.getDataSet("distances");
    distance_data_set.write(distance_matrix);
    std::cout << "Wrote distances to: " << arguments.out_path.string() << std::endl;
  }

  return EXIT_SUCCESS;
}