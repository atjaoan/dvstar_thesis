#pragma once

#include <filesystem>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

#include "cluster_container.hpp"
#include "vlmc_container.hpp"
#include "distances/dvstar.hpp"
#include "global_aliases.hpp"

namespace parser {

enum Mode {
  compare,
  kmer_major
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

enum VLMC_Rep {
  vlmc_vector, 
  vlmc_sorted_vector,
  vlmc_b_tree,
  vlmc_hashmap,
  vlmc_veb,
  vlmc_ey,
  vlmc_alt_btree,
  vlmc_sorted_search,
  vlmc_sorted_search_ey
};

struct cli_arguments {
  Mode mode{Mode::compare};
  Distance_function dist_fn{Distance_function::dvstar};
  std::filesystem::path first_VLMC_path{};
  std::filesystem::path second_VLMC_path{};
  std::filesystem::path out_path{};
  size_t dop {1};
  int set_size {-1};
  Cluster_Rep cluster{Cluster_Rep::cluster_vector};
  VLMC_Rep vlmc{VLMC_Rep::vlmc_vector}; 
  size_t background_order {0}; 
};

template <typename VC>
std::function<out_t(VC &, VC &)> parse_distance_function(cli_arguments arguments) {
  if (arguments.dist_fn == parser::Distance_function::dvstar) {
    auto fun = [&](auto &left, auto &right) {
      return distance::dvstar<VC>(left, right, arguments.background_order);
    };
    return fun; 
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
      {"compare", Mode::compare},
      {"kmer-major", Mode::kmer_major},
  };

  std::map<std::string, Distance_function> function_map{
      {"kl", Distance_function::kl},
      {"dvstar", Distance_function::dvstar},
  };

  std::map<std::string, VLMC_Rep> VLMC_Rep_map{
      {"vector", VLMC_Rep::vlmc_vector},
      {"sorted-vector", VLMC_Rep::vlmc_sorted_vector},
      {"b-tree", VLMC_Rep::vlmc_b_tree},
      {"hashmap", VLMC_Rep::vlmc_hashmap},
      {"veb", VLMC_Rep::vlmc_veb},
      {"ey", VLMC_Rep::vlmc_ey},
      {"alt-btree", VLMC_Rep::vlmc_alt_btree},
      {"sorted-search", VLMC_Rep::vlmc_sorted_search},
      {"sorted-search-ey", VLMC_Rep::vlmc_sorted_search_ey}
  };

  std::map<std::string, Cluster_Rep> cluster_rep_map{
      {"vector", Cluster_Rep::cluster_vector}
  };

  app.add_option(
         "-m,--mode", arguments.mode,
         "Program mode, 'compare', 'kmer-major'." "Place holder descriptiopn")
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
                 "Path to hdf5 file where scores will be stored.");
  
  app.add_option("-n,--max-dop", arguments.dop,
                 "Degree of parallelism. Default 1 (sequential).");

  app.add_option("-v,--vlmc-rep", arguments.vlmc,
                 "Vlmc container representation to use.")
      ->transform(CLI::CheckedTransformer(VLMC_Rep_map, CLI::ignore_case));

  app.add_option("-c,--cluster-rep", arguments.cluster,
                 "Cluster container representation to use.")
      ->transform(CLI::CheckedTransformer(cluster_rep_map, CLI::ignore_case));

  app.add_option("-b,--background-order", arguments.background_order,
                 "Background order.");

  app.add_option("-a, --set-size", arguments.set_size,
                    "Number of VLMCs to compute distance function on.");
}
}
