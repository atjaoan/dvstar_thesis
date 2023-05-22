# Comparing variable-length Markov chains
This repository contains the code for comparing a directory/directories of VLMCs. It builds upon the existing repositories [dvstar](https://github.com/Schlieplab/dvstar) and [pst-classifier-seqan repository](https://github.com/Schlieplab/PstClassifierSeqan) by Joel Gustafsson. The VLMCs comparisons are speedup by a factor of 2-25 compared to the pst-classifier-seqan repository. The current implementation requires that the VLMCs are built using Joels implementation since the VLMCs are expected to be of binary format and with a certain structure. 

## Publication

Publication pending...

## Compilation

To compile the program, please first download all submodules:

```shell script
git submodule update --init --recursive
```

Create a build directory, configure with cmake and build:

```shell script
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=Release ..
make
```

## Execution

This provides an executable `dist`, which can be used as follows:

```shell
% ./dist --help
Comparisons of variable-length Markov chains.
Usage: ./dist [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -p,--VLMC-path TEXT         Required for distance calculation. 'Primary' path to saved bintree directory. If '-s' is empty it will compute the inter-distance between the trees of the directory.
  -s,--snd-VLMC-path TEXT     Optional 'Secondary' path to saved bintree directory. Calculates distance between the trees specified in -p (primary) and -s (secondary).
  -o,--matrix-path TEXT       Path to hdf5 file where scores will be stored. If left empty, distances will be printed to shell.
  -n,--max-dop UINT           Degree of parallelism. Default 1 (sequential).
  -v,--vlmc-rep               VLMC container to use for comparison, see paper for more details. If unsure use standard (sorted-search). Available options: 'sorted-search', 'sorted-vector', 'b-tree', 'ey', 'hashmap', 'kmer-major', 'veb', 'vector'
                              Vlmc container representation to use.
  -b,--background-order UINT  Background order.
  -a,--set-size INT           Number of VLMCs to compute distance function on. If left empty will load all VLMCs in the 'primary' and 'secondary' directories. Otherwise, loads the specified amount from the given directories. 
```

For example, to compare two directories of VLMCs using 8 cores, run (from build/):

```shell
./dist --VLMC-path ../tests/dir_p --snd-VLMC-path ../tests/dir_s --max-dop 8
```

## Headers

If, for some reason, you wanted to include the code in some other project, this directory can be included with CMAKE as
`add_subdirectory(/path/to/dvstar_thesis)`, or alternatively, include the headers
in `/path/to/dvstar_thesis/include`.
Note that with the second option, the submodules need to be included manually.

### Container solution / apptainer (previously singularity)

We provide an [apptainer container](https://apptainer.org/). The definition file needs to be modified to
include the path to kmc, see the fifth and sixth line. 

The submodule CLI11 needs to be adjusted as well, specifically the folder `tests/mesonTest/subprojects/CLI11` should be removed. 

The container is built by running:

```shell script
apptainer build --fakeroot thesis.sif thesis.def
```

Using the created image, programs can be run by either using the ```apptainer exec ./program``` or ```apptainer run``` to run the programs in ```%runscripts``` of the .def file. It is also possible to spawn a shell in the container with 
```shell script
apptainer shell thesis.sif
```

To output benchmarking results, the ```--bind``` option must be specified such that the container has the output directory mounted (the containers own file system is read-only). First, create two empty directory named ```csv_results``` and ```hdf5_results``` in the same directory as the .sif file is located. 
```shell script
mkdir csv_results hdf5_results
```
Benchmarking then is done with
```shell script
apptainer run --bind ./csv_results:/thesis/csv_results/,./hdf5_results:/thesis/hdf5_results thesis.sif
```