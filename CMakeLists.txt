cmake_minimum_required (VERSION 3.4...3.12)

project (kmc_vlmc CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Threads REQUIRED)

include_directories("${CMAKE_SOURCE_DIR}/KMC/kmc_api")

add_library(kmc_file SHARED KMC/kmc_api/kmc_file.cpp)
add_library(kmer_api SHARED KMC/kmc_api/kmer_api.cpp)
add_library(mmer SHARED KMC/kmc_api/mmer.cpp)
set(KMC_Libraries TRUE)

# STXXL
add_subdirectory(stxxl)
# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})


add_executable(build_vlmc kmc_vlmc/src/build_vlmc.cpp)
target_link_libraries(build_vlmc Threads::Threads)
target_link_libraries(build_vlmc kmc_file kmer_api mmer)
target_link_libraries(build_vlmc ${STXXL_LIBRARIES})


add_subdirectory(kmc_vlmc/tests)