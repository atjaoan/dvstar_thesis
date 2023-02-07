#include <filesystem>

#include "vlmc_from_kmers/context_archive.hpp"
#include "vlmc_from_kmers/kmer_container.hpp"
#include "vlmc_from_kmers/similarity_pruning.hpp"
#include "vlmc_from_kmers/kmer.hpp"

constexpr int max_k = 255;

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;
//KmerContainer<class Comparator = ReverseKMerComparator<max_k>>;

namespace get_trees{
  
bool has_next(std::ifstream &file_stream) { return file_stream.peek() != EOF; }

vlmc::VLMCKmer &next(cereal::BinaryInputArchive &iarchive, vlmc::VLMCKmer &kmer) {
  iarchive(kmer);
  return kmer;
}

void get_trees(const std::filesystem::path &directory){
  // vlmc::HashMapKmerContainer kmer_con(); 
  std::vector<std::vector<vlmc::VLMCKmer>> trees{};

  for (const auto& dir_entry : recursive_directory_iterator(directory)) {
    std::ifstream ifs(dir_entry.path(), std::ios::binary);
    cereal::BinaryInputArchive archive(ifs);

    std::vector<vlmc::VLMCKmer> tree{};
    vlmc::VLMCKmer kmer{};

    while (has_next(ifs)){
      next(archive, kmer);
      tree.push_back(kmer);
    }
    trees.push_back(tree);
    ifs.close();
  } 

  std::cout << trees.size() << std::endl; 
  
}
}
