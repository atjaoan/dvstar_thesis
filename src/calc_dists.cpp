#include <Eigen/Dense>

#include "optimize_cmp_kmers/parser.hpp"
#include "optimize_cmp_kmers/get_trees.hpp"
#include "optimize_cmp_kmers/cluster_container.hpp"
#include "optimize_cmp_kmers/vlmc_template.hpp"

// Distance Functions 
#include "optimize_cmp_kmers/distances/dvstar.hpp"
#include "optimize_cmp_kmers/distances/kl_divergence.hpp"

using matrix_t = Eigen::MatrixXd;

using vlmc_c = container::VLMC_vector;
using cluster_c = cluster::Cluster_vector<vlmc_c>;

std::function<float(vlmc_c &, vlmc_c &)>
parse_distance_function(parser::Distance_function dist_fn) {

  if (dist_fn == parser::Distance_function::dvstar) {
    return distance::dvstar<vlmc_c>; 
  } 
  else if (dist_fn ==  parser::Distance_function::kl) {
    return distance::kl<vlmc_c>; 
  }  
  throw std::invalid_argument("Invalid distance function name.");
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
  if(arguments.mode == parser::Mode::compare){
    if(arguments.in_path.empty()){
      std::cerr
          << "Error: A input path to .bintree files has to be given for comparison operation."
          << std::endl;
      return EXIT_FAILURE;
    }

    if(arguments.to_path.empty()){
      cluster_c trees{}; 
      get_trees::get_trees<cluster_c, vlmc_c>(arguments.in_path, trees);

    } else {
      cluster_c left_trees{};
      cluster_c right_trees{};
      get_trees::get_trees<cluster_c, vlmc_c>(arguments.in_path, left_trees);
      get_trees::get_trees<cluster_c, vlmc_c>(arguments.to_path, right_trees);
    }
  }

  return EXIT_SUCCESS;
}