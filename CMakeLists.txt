cmake_minimum_required (VERSION 3.4...3.12)

project (CountVLMC CXX)

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

include_directories("${CMAKE_SOURCE_DIR}/KMC/kmc_api")

add_library(kmc_file STATIC KMC/kmc_api/kmc_file.cpp)
add_library(kmer_api STATIC KMC/kmc_api/kmer_api.cpp)
add_library(mmer STATIC KMC/kmc_api/mmer.cpp)
set(KMC_Libraries TRUE)

# STXXL
set(USE_BOOST OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/stxxl)
# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})

# Parallelization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbb -fopenmp")

# CLI11 for command line interface
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLI11)

set(SKIP_PERFORMANCE_COMPARISON ON)

# Cereal for serialisation
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cereal)

set(CountVLMC_LIBRARIES ${STXXL_LIBRARIES} tbb kmc_file kmer_api mmer CLI11::CLI11 cereal::cereal)

add_library(CountVLMC INTERFACE)
add_library(CountVLMC::CountVLMC ALIAS CountVLMC)

target_link_libraries(CountVLMC INTERFACE ${CountVLMC_LIBRARIES})
target_include_directories(CountVLMC INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/kmc_vlmc/src/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KMC/kmc_api/>
        $<INSTALL_INTERFACE:kmc_vlmc/src/>)

if(MAIN_PROJECT)
    add_executable(build_vlmc kmc_vlmc/src/build_vlmc.cpp)
    target_link_libraries(build_vlmc Threads::Threads)
    target_link_libraries(build_vlmc ${CountVLMC_LIBRARIES})
endif()



#add_subdirectory(kmc_vlmc/tests)
