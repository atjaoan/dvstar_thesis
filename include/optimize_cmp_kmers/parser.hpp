#pragma once

#include <filesystem>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

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
  std::filesystem::path first_VLMC_path;
  std::filesystem::path second_VLMC_path;
  std::filesystem::path out_path{};
};

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
}



}
