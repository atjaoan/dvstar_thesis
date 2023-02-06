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
  std::filesystem::path in_path;
  std::filesystem::path to_path;
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
      "--in-path", arguments.in_path,
      "Path to saved bintree directory.");

  app.add_option(
      "--to-path", arguments.to_path,
      "Path to saved bintree directory.");

  app.add_option("-o,--out-path", arguments.out_path,
                 "Path to matrix of distances.");
}



}
