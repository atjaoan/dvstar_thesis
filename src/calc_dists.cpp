#pragma once 

#include <highfive/H5File.hpp>

#include "parser.hpp"
#include "get_cluster.hpp"
#include "calc_dists.hpp"
#include "global_aliases.hpp"
#include "utils.hpp"

template <typename VC>
matrix_t calculate_cluster_distance(parser::cli_arguments arguments, const size_t nr_cores){
  auto distance_function = parser::parse_distance_function<VC>(arguments);
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
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_sorted_vector){
    return calculate_cluster_distance<container::VLMC_sorted_vector>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_b_tree){
    return calculate_cluster_distance<container::VLMC_B_tree>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_hashmap){
    return calculate_cluster_distance<container::VLMC_hashmap>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_veb){
    return calculate_cluster_distance<container::VLMC_Veb>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_ey){
    return calculate_cluster_distance<container::VLMC_Eytzinger>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_alt_btree){
    return calculate_cluster_distance<container::VLMC_Alt_Btree>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_sorted_search){
    return calculate_cluster_distance<container::VLMC_sorted_search>(arguments, nr_cores);
  } else if (vlmc_container==parser::VLMC_Rep::vlmc_sorted_search_ey){
    return calculate_cluster_distance<container::VLMC_sorted_search_ey>(arguments, nr_cores);
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

std::map<parser::VLMC_Rep, std::string> VLMC_Rep_map{
      {parser::VLMC_Rep::vlmc_vector, "vector"},
      {parser::VLMC_Rep::vlmc_sorted_vector, "sorted-vector"},
      {parser::VLMC_Rep::vlmc_b_tree, "b-tree"},
      {parser::VLMC_Rep::vlmc_hashmap, "hashmap"},
      {parser::VLMC_Rep::vlmc_veb, "veb"}, 
      {parser::VLMC_Rep::vlmc_ey, "ey"},
      {parser::VLMC_Rep::vlmc_alt_btree, "alt-btree"},
      {parser::VLMC_Rep::vlmc_sorted_search, "sorted-search"},
      {parser::VLMC_Rep::vlmc_sorted_search_ey, "sorted-search-ey"}
  };

std::string get_group_name(parser::cli_arguments arguments){
  auto ret_str = "dop-" + std::to_string(arguments.dop); 
  ret_str += "-set-size-" + std::to_string(arguments.set_size); 
  ret_str += "-bo-" + std::to_string(arguments.background_order);
  if (arguments.mode == parser::Mode::compare) {
    auto in_data = utils::get_filename(arguments.first_VLMC_path);
    return VLMC_Rep_map[arguments.vlmc] + "-" + in_data + "-" + ret_str; 
  } else {
    auto in_data = utils::get_filename(arguments.first_VLMC_path);
    return "kmer-major-" + in_data + "-" + ret_str;
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

  matrix_t distance_matrix;

  if(arguments.mode==parser::Mode::compare){
    std::cout << "Running Normal implementation" << std::endl; 
    distance_matrix = apply_container(arguments, arguments.vlmc, nr_cores);
  } else {
    std::cout << "Running Kmer-Major implementation" << std::endl; 
    distance_matrix = calculate_kmer_major(arguments, nr_cores);
  }

  if (arguments.out_path.empty()) {
    utils::print_matrix(distance_matrix);
  } 
  else if (arguments.out_path.extension() == ".h5" ||
             arguments.out_path.extension() == ".hdf5") {
    HighFive::File file{arguments.out_path, HighFive::File::OpenOrCreate};
    auto group_name = get_group_name(arguments); 

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