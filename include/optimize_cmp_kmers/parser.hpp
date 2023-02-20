#pragma once

#include <filesystem>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

#include "cluster_container.hpp"
#include "vlmc_container.hpp"

// Distance Functions 
#include "distances/dvstar.hpp"
#include "distances/kl_divergence.hpp"

namespace parser {

using cluster_c = container::Cluster_Container;
using vlmc_c = container::VLMC_Container;

enum Mode {
  compare
};

enum Distance_function {
  d2, 
  d2star, 
  dvstar, 
  nearest_dvstar, 
  penalized_dvstar, 
  kl, 
  kl_both, 
  nll, 
  nll_background,
  cv,
  cv_estimation
};

enum Cluster_Rep {
  cluster_vector
};

enum Vlmc_Rep {
  vlmc_vector, 
  vlmc_multi_vector,
  vlmc_sorted_vector,
  vlmc_b_tree,
  vlmc_hashmap
};

struct cli_arguments {
  Mode mode{Mode::compare};
  Distance_function dist_fn{Distance_function::dvstar};
  std::filesystem::path first_VLMC_path{};
  std::filesystem::path second_VLMC_path{};
  std::filesystem::path out_path{};
  size_t dop {1};
  Cluster_Rep cluster{Cluster_Rep::cluster_vector};
  Vlmc_Rep vlmc{Vlmc_Rep::vlmc_vector}; 
  size_t background_order {0}; 
};

std::function<double(vlmc_c &, vlmc_c &)>
parse_distance_function(parser::Distance_function dist_fn, size_t background_order) {
  if (dist_fn == parser::Distance_function::dvstar) {
    auto fun = [&](auto &left, auto &right) {
      return distance::dvstar(left, right, background_order);
    };
    return fun; 
  } 
  else if (dist_fn ==  parser::Distance_function::kl) {
    return distance::kl; 
  }  
  throw std::invalid_argument("Invalid distance function name.");
}

size_t parse_dop(size_t requested_cores){
  if(requested_cores < 1){
    throw std::invalid_argument("Too low degree of parallelism, must be >= 1");
  } else{
    return requested_cores;
  }
}

void add_options(CLI::App &app, cli_arguments &arguments) {
  std::map<std::string, Mode> mode_map{
      {"compare", Mode::compare}};

  std::map<std::string, Distance_function> function_map{
      {"kl", Distance_function::kl},
      {"dvstar", Distance_function::dvstar},
  };

  std::map<std::string, Vlmc_Rep> vlmc_rep_map{
      {"vector", Vlmc_Rep::vlmc_vector},
      {"multi-vector", Vlmc_Rep::vlmc_multi_vector},
      {"sorted-vector", Vlmc_Rep::vlmc_sorted_vector},
      {"b-tree", Vlmc_Rep::vlmc_b_tree},
      {"hashmap", Vlmc_Rep::vlmc_hashmap},
  };

  std::map<std::string, Cluster_Rep> cluster_rep_map{
      {"vector", Cluster_Rep::cluster_vector}
  };

  app.add_option(
         "-m,--mode", arguments.mode,
         "Program mode, 'compare'." "Place holder descriptiopn")
      ->transform(CLI::CheckedTransformer(mode_map, CLI::ignore_case));

  app.add_option("--function", arguments.dist_fn,
                 "Distance function to use 'dvstar', or 'kl'.")
      ->transform(CLI::CheckedTransformer(function_map, CLI::ignore_case));

  app.add_option(
      "-p,--VLMC-path", arguments.first_VLMC_path,
      "Required for distance calcualtion. Path to saved bintree directory. Computes the inter-distance between the trees of the directory.");

  app.add_option(
      "-s,--snd-VLMC-path", arguments.second_VLMC_path,
      "Optional path to saved bintree directory. Calculates distance between the trees in VLMC-path and snd-VLMC-path.");

  app.add_option("-o,--matrix-path", arguments.out_path,
                 "Path to matrix of distances.");
  
  app.add_option("-n,--max-dop", arguments.dop,
                 "Degree of parallelism. Default 1 (sequential).");

  app.add_option("-v,--vlmc-rep", arguments.vlmc,
                 "Vlmc container representation to use.")
      ->transform(CLI::CheckedTransformer(vlmc_rep_map, CLI::ignore_case));

  app.add_option("-c,--cluster-rep", arguments.cluster,
                 "Cluster container representation to use.")
      ->transform(CLI::CheckedTransformer(cluster_rep_map, CLI::ignore_case));

  app.add_option("-b,--background-order", arguments.background_order,
                 "Background order.");
}

}
