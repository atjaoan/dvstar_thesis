cmake_minimum_required (VERSION 3.4...3.12)

project (CountVLMC CXX)

add_definitions(-w)

option(CountVLMC_COMPILE_TESTS "To compile tests, specify -DCountVLMC_COMPILE_TESTS=ON")
option(COMPILE_OLD "To compile original dvstar implementation, specify -DCOMPILE_OLD=ON")

set(MAIN_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(MAIN_PROJECT ON)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "/usr/share/cmake/Modules/")

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED CXX)
find_package(TBB REQUIRED)
find_package(unordered_dense CONFIG REQUIRED)

add_library(kmc_file STATIC submodules/KMC/kmc_api/kmc_file.cpp)
add_library(kmer_api STATIC submodules/KMC/kmc_api/kmer_api.cpp)
add_library(mmer STATIC submodules/KMC/kmc_api/mmer.cpp)
set(KMC_Libraries TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
# STXXL
set(USE_BOOST OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/submodules/stxxl)
# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# Parallelization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbb -fopenmp -fpermissive")

# CLI11 for command line interface
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/submodules/CLI11)

set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_USE_EIGEN ON)
set(HIGHFIVE_EXAMPLES FALSE)

# HighFive
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/submodules/HighFive)

set(SKIP_PERFORMANCE_COMPARISON ON)

# Add include paths
set(INCLUDE_DIRS ${STXXL_INCLUDE_DIRS} "${PROJECT_SOURCE_DIR}/submodules/stxxl/include" "${PROJECT_SOURCE_DIR}/submodules/robin-hood-hashing/src/include" "${PROJECT_SOURCE_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}/submodules/KMC/kmc_api"
                                         "${CMAKE_CURRENT_SOURCE_DIR}/submodules/cereal/include" "${CMAKE_CURRENT_SOURCE_DIR}/submodules/eigen" "${CMAKE_CURRENT_SOURCE_DIR}/include/optimize_cmp_kmers"
                                         "${PROJECT_SOURCE_DIR}/submodules/unordered_dense/include/ankerl" "${CMAKE_CURRENT_SOURCE_DIR}/submodules/thread-pool")

include_directories("${INCLUDE_DIRS}")

if(NOT MAIN_PROJECT)
    set(CountVLMC_INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE)
endif()

set(CountVLMC_LIBRARIES ${STXXL_LIBRARIES} tbb kmc_file kmer_api mmer CLI11::CLI11 unordered_dense::unordered_dense HighFive)

add_library(CountVLMC INTERFACE)
add_library(CountVLMC::CountVLMC ALIAS CountVLMC)

target_link_libraries(CountVLMC INTERFACE ${CountVLMC_LIBRARIES})
target_include_directories(CountVLMC INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KMC/kmc_api/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>
        $<INSTALL_INTERFACE:include/>)

if(MAIN_PROJECT)
    if(COMPILE_OLD)
        add_executable(dvstar src/build_vlmc.cpp include/vlmc_from_kmers/context_archive.hpp include/vlmc_from_kmers/estimators.hpp include/vlmc_from_kmers/sequencing_adjustment.hpp)
        target_link_libraries(dvstar Threads::Threads)
        target_link_libraries(dvstar ${CountVLMC_LIBRARIES})
    endif()
    add_executable(dist src/calc_dists.cpp)
    add_executable(time_benchmark tests/our_tests/time_benchmark.cpp)
    target_link_libraries(dist ${CountVLMC_LIBRARIES})
    target_link_libraries(time_benchmark ${CountVLMC_LIBRARIES})

    if (CountVLMC_COMPILE_TESTS)
        add_subdirectory(tests)
    endif()
endif()
