#include "optimize_cmp_kmers/parser.hpp"

#include "optimize_cmp_kmers/distances/dvstar.hpp"
#include "optimize_cmp_kmers/distances/kl_divergence.hpp"

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
    //vlmc::dump_path(arguments.in_path, arguments.out_path);
  }

  std::cout << arguments.dist_fn << std::endl;
  return EXIT_SUCCESS;
}