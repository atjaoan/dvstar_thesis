#pragma once 

#include <Eigen/Core>
#include "vlmc_from_kmers/kmer.hpp"

// OUTPUT TYPE
using out_t = double;

// EIGEN 
using eigen_t = Eigen::Array4d;
using eigenx_t = Eigen::ArrayX4d;  
using matrix_t = Eigen::MatrixXd;

// KMER 
using Kmer = vlmc::VLMCKmer;

// FILESYSTEM 
using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;

// ERROR TOLERANCE
const out_t error_tolerance = 1E-8;