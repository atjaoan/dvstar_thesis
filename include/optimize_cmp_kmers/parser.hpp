#pragma once

#include <filesystem>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

// Distance Functions 
#include "optimize_cmp_kmers/distances/dvstar.hpp"
#include "optimize_cmp_kmers/distances/kl_divergence.hpp"

namespace parser {

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
  std::filesystem::path in_path;
  std::filesystem::path to_path;
  std::filesystem::path out_path{};
};

template <typename VLMC_Container>
std::function<float(VLMC_Container &, VLMC_Container &)>
parse_distance_function(parser::Distance_function dist_fn) {

  if (dist_fn == parser::Distance_function::dvstar) {
    return distance::dvstar<VLMC_Container>; 
  } 
  else if (dist_fn ==  parser::Distance_function::kl) {
    return distance::kl<VLMC_Container>; 
  }  
  throw std::invalid_argument("Invalid distance function name.");
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
      "--in-path", arguments.in_path,
      "Path to saved bintree directory.");

  app.add_option(
      "--to-path", arguments.to_path,
      "Path to saved bintree directory.");

  app.add_option("-o,--out-path", arguments.out_path,
                 "Path to matrix of distances.");
}



}
