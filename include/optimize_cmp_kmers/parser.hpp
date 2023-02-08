#pragma once

#include <filesystem>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

#include "vlmc_container.hpp"

// Distance Functions 
#include "distances/dvstar.hpp"
#include "distances/kl_divergence.hpp"

namespace parser {

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

struct cli_arguments {
  Mode mode{Mode::compare};
  Distance_function dist_fn{Distance_function::dvstar};
  std::filesystem::path first_VLMC_path{};
  std::filesystem::path second_VLMC_path{};
  std::filesystem::path out_path{};
  size_t dop {1};
};

std::function<float(vlmc_c &, vlmc_c &)>
parse_distance_function(parser::Distance_function dist_fn) {

  if (dist_fn == parser::Distance_function::dvstar) {
    return distance::dvstar; 
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
}

}
