# Building a variable-length Markov chain from k-mers
This repository contains the code for constructing a variable-length Markov chain (VLMC) from a lexicographically
sorted list of _k_-mers.  This enables VLMCs to be constructed on even very large genomes, as well as collections
of sequences.

This repository deals only with the task of constructing the VLMC, for scoring of sequences, or similarities between
two VLMCs, please see the [pst-classifier-seqan repository](https://github.com/Schlieplab/PstClassifierSeqan).

## Publication
Publication pending...

## Compilation
To compile the program, please first download all submodules:

```shell script
git submodule update --init --recursive
```

### Container solution / Singularity
We provide a [singularity container](https://sylabs.io/singularity).  The definition file may need to be modified to 
include the path to kmc, see the fifth and sixth line. The container is built by running:   

```shell script
sudo singularity build vlmc-from-kmers.sif vlmc-from-kmers.def
```

### Manually
Create a build directory, configure with cmake and build:

```shell script
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=Release ..
make
```

__Note that kmc3 needs to be installed separately, and be in the current working directory, or in your path.__ Both `kmc` and `kmc_tools` are required. 

This provides an executable `build_vlmc`, which can be used as follows:

```shell
./build_vlmc --help
Variable-length Markov chain construction construction using k-mer counter.
Usage: ./build_vlmc [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -m,--mode ENUM:value in {build->0,dump->2,score->1,bic->3} OR {0,2,1,3}
                              Program mode, 'build', 'dump', or 'score'.
  -p,--fasta-path TEXT        Path to fasta file.  Required for 'build' and 'score' modes.
  --in-path TEXT              Path to saved tree file.  Required for 'dump' and 'score' modes.
  -o,--out-path TEXT          Path to output file.  The VLMCs are stored as binary, and can be read by the 'dump' or 'score' modes.  Required for 'build' and 'dump' modes.
  -t,--temp-path TEXT         Path to temporary folder for the external memory algorithms.  For good performance, this needs to be on a local machine.  For sorting, at least 2GB will be allocated to this path.  Defaults to ./tmp
  -c,--min-count INT          Minimum count required for every k-mer in the tree.
  -k,--threshold FLOAT        Kullback-Leibler threshold.
  -d,--max-depth INT          Maximum depth for included k-mers.
  -a,--pseudo-count-amount FLOAT
                              Size of pseudo count for probability estimation. See e.g. https://en.wikipedia.org/wiki/Additive_smoothing .
  -i,--in-or-out-of-core ENUM:value in {external->0,hash->2,internal->1} OR {0,2,1}
                              Specify 'internal' for in-core or 'external for out-of-core memory model.  Out of core is slower, but is not memory bound.

```

For example, to construct a VLMC, run:

```shell
./build_vlmc --fasta-path NC_022098.1.fasta --threshold 3.9075 --max-depth 15 --min-count 100 --out-path NC_022098.1.bintree --temp-path tmp
```

To view the contents of the VLMC, run:
```shell
./build_vlmc --mode dump --in-path NC_022098.1.bintree
```

## Headers
If, for some reason, you wanted to include the code in some other project, this directory can be included with CMAKE as 
`add_subdirectory(/path/to/vlmc-from-kmers)`, or alternatively, include the headers in `/path/to/vlmc-from-kmers/include`.
Note that with the second option, the submodules need to be included manually.
