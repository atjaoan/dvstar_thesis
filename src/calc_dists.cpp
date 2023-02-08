#include "parser.hpp"
#include "get_trees.hpp"
#include "cluster_container.hpp"
#include "vlmc_template.hpp"
#include "calc_dists.hpp"

using matrix_t = Eigen::MatrixXd;

using vlmc_c = container::VLMC_vector;
using cluster_c = cluster::Cluster_vector<vlmc_c>;


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
    if(arguments.first_VLMC_path.empty()){
      std::cerr
          << "Error: A input path to .bintree files has to be given for comparison operation."
          << std::endl;
      return EXIT_FAILURE;
    }
    auto distance_function = parser::parse_distance_function<vlmc_c>(arguments.dist_fn);

    if(arguments.second_VLMC_path.empty()){
      cluster_c trees{};
      get_trees::get_trees<cluster_c, vlmc_c>(arguments.first_VLMC_path, trees);
      matrix_t distance_matrix;
    } else {
      cluster_c left_trees{};
      cluster_c right_trees{};
      get_trees::get_trees<cluster_c, vlmc_c>(arguments.first_VLMC_path, left_trees);
      get_trees::get_trees<cluster_c, vlmc_c>(arguments.second_VLMC_path, right_trees);
      matrix_t distance_matrix;
    }
  }

  return EXIT_SUCCESS;
}