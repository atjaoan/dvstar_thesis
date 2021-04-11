#pragma once

std::tuple<std::string, size_t, std::array<size_t, 4>>
read_line(std::string &line) {
  std::istringstream input(line);
  std::vector<std::string> values;
  for (std::string c; std::getline(input, c, ' ');) {
    values.push_back(c);
  }

  return {values[0], std::stoull(values[1]),
          std::array<size_t, 4>{std::stoull(values[2]), std::stoull(values[3]),
                                std::stoull(values[4]),
                                std::stoull(values[5])}};
}

std::tuple<std::string, size_t, std::array<size_t, 4>>
read_tree_line(std::string &line) {
  std::istringstream input(line);
  std::vector<std::string> values;
  for (std::string c; std::getline(input, c, ' ');) {
    values.push_back(c);
  }

  std::string kmer = values[2];
  if (kmer == "#") {
    kmer = "";
  }

  return {kmer, std::stoull(values[9]),
          std::array<size_t, 4>{
              std::stoull(values[11]) - 1, std::stoull(values[12]) - 1,
              std::stoull(values[13]) - 1, std::stoull(values[14]) - 1}};
}