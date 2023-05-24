cmake_minimum_required (VERSION 3.4...3.12)

project (CountVLMC CXX)

add_definitions(-w)

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

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")

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
set(INCLUDE_DIRS ${STXXL_INCLUDE_DIRS} "${PROJECT_SOURCE_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}/submodules/cereal/include" "${CMAKE_CURRENT_SOURCE_DIR}/submodules/eigen" 
                                        "${PROJECT_SOURCE_DIR}/submodules/unordered_dense/include/ankerl" "${CMAKE_CURRENT_SOURCE_DIR}/submodules/thread-pool")

include_directories("${INCLUDE_DIRS}")

if(NOT MAIN_PROJECT)
    set(CountVLMC_INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE)
endif()

set(CountVLMC_LIBRARIES tbb CLI11::CLI11 unordered_dense::unordered_dense HighFive)

add_library(CountVLMC INTERFACE)
add_library(CountVLMC::CountVLMC ALIAS CountVLMC)

target_link_libraries(CountVLMC INTERFACE ${CountVLMC_LIBRARIES})
target_include_directories(CountVLMC INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>
        $<INSTALL_INTERFACE:include/>)

if(MAIN_PROJECT)
    add_executable(dist src/calc_dists.cpp)
    target_link_libraries(dist ${CountVLMC_LIBRARIES})
endif()
