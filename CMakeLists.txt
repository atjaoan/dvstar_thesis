cmake_minimum_required (VERSION 3.4...3.12)

project (kmc_vlmc CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "/usr/share/cmake/Modules/")

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED CXX)
find_package(TBB REQUIRED)

include_directories("${CMAKE_SOURCE_DIR}/KMC/kmc_api")

add_library(kmc_file STATIC KMC/kmc_api/kmc_file.cpp)
add_library(kmer_api STATIC KMC/kmc_api/kmer_api.cpp)
add_library(mmer STATIC KMC/kmc_api/mmer.cpp)
set(KMC_Libraries TRUE)

# STXXL
add_subdirectory(stxxl)
# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})

# Parallelization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbb -fopenmp")

# CLI11 for command line interface
add_subdirectory(CLI11)

# Cereal for serialisation
include_directories(cereal/include)

add_executable(build_vlmc kmc_vlmc/src/build_vlmc.cpp)
target_link_libraries(build_vlmc Threads::Threads)
target_link_libraries(build_vlmc kmc_file kmer_api mmer)
target_link_libraries(build_vlmc ${STXXL_LIBRARIES})
target_link_libraries(build_vlmc CLI11::CLI11)
target_link_libraries(build_vlmc tbb)


add_subdirectory(kmc_vlmc/tests)
