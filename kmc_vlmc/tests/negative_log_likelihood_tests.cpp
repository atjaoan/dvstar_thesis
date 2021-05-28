#include <gtest/gtest.h>

#include "../src/build_vlmc.hpp"
#include "../src/negative_log_likelihood.hpp."

#include "read_helper.hpp"

using namespace vlmc;

class NegativeLogLikelihoodTests : public ::testing::Test {
protected:
  void SetUp() override {}
};

void generate_fasta(const std::filesystem::path &path) {
  std::ofstream file_stream{path};
  file_stream << ">Test" << std::endl;
  file_stream << "ACGTACGATCGTACGATACGTACGATCGTACGATACGTACGATCGT" << std::endl;
  file_stream << "CCGTACGATCGTACGAGACGTACGATCGTACGATCCGTACTATCGT" << std::endl;
  file_stream << "TTCAGCTAGTCAGATCTAGCTAGCTAGCTACGTACTAGCTAGCGAT" << std::endl;
  file_stream << "ACGATCGTAGTCGATCAGTCATCTACTGACACGAGCGCTCGATGAC" << std::endl;
  file_stream << "ATCATCATCATAATTCAGTATATAGCATTATCGATTATCGATTACG" << std::endl;
  file_stream << "ATCGATCGGATCGATCGATCGATCGATCGATCGATCGATCGATCGA" << std::endl;
  file_stream.close();
}

TEST_F(NegativeLogLikelihoodTests, Correct) {
  std::filesystem::path temp_fasta_1{"tmp1.fasta"};
  std::filesystem::path temp_fasta_2{"tmp2.fasta"};
  std::filesystem::path temp{"tmp"};
  std::filesystem::path out{"_test.tree"};

  generate_fasta(temp_fasta_1);
  generate_fasta(temp_fasta_2);

  build(temp_fasta_1, 2, 1, 3.9075, out, temp, "internal");

  double nll = negative_log_likelihood(temp_fasta_2, temp, out, "internal", 2);

  EXPECT_EQ(0.0, nll);
}