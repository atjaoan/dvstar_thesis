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

### Container solution / apptainer (previously singularity)
We provide a [apptainer container](https://apptainer.org/).  The definition file may need to be modified to
include the path to kmc, see the fifth and sixth line. The container is built by running:   

```shell script
sudo apptainer build vlmc-from-kmers.sif vlmc-from-kmers.def
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
% ./build_vlmc --help
Variable-length Markov chain construction construction using k-mer counter.
Usage: ./build_vlmc [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -m,--mode ENUM:value in {bic->3,build->0,build-from-kmc-db->4,dump->2,dvstar->5,score->1} OR {3,0,4,2,5,1}
                              Program mode, 'build', 'build-from-kmc-db', 'dump', 'score', or 'dvstar'.  For build-from-kmc-db, the kmc db needs to include all k-mers (not in canonical form), with no minimum count cutoff.  The length of the k-mers needs to be set to 1 more than the maximum depth of the VLMC.  The kmc db also has to be sorted.
  -p,--fasta-path TEXT        Path to fasta file.  Required for 'build' and 'score' modes.
  --in-path TEXT              Path to saved tree file or kmc db file.  Required for 'build-from-kmc-db', 'dump', 'score', and 'dvstar' modes.  For 'build-from-kmc-db', the kmc db file needs to be supplied without the file extension.
  --to-path TEXT              Path to saved tree file.  Required for 'dvstar' mode.
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
./build_vlmc --fasta-path NC_022098.1.fasta --threshold 3.9075 --max-depth 4 --min-count 100 --out-path NC_022098.1.bintree --temp-path tmp
```

To view the contents of the VLMC, run:
```shell
./build_vlmc --mode dump --in-path NC_022098.1.bintree
```

To compute the dvstar similarity between two VLMCs:
```shell
./build_vlmc --mode dvstar --in-path NC_022098.1.bintree --to-path NC_022098.1.bintree
```

## Visualisation

To produce a graphical representation of the VLMCs, see the [`visualiser.py`](visualiser.py) script.
It needs the `plotly` and `typer` python packages to run, which can be installed e.g. through `pip install --user plotly typer`.

```shell
python visualiser.py NC_022098.1.bintree
```

Produces an interactive visualisation (should open in your browser) as well as a pdf and png file in the `images` directory.

## Headers
If, for some reason, you wanted to include the code in some other project, this directory can be included with CMAKE as 
`add_subdirectory(/path/to/vlmc-from-kmers)`, or alternatively, include the headers in `/path/to/vlmc-from-kmers/include`.
Note that with the second option, the submodules need to be included manually.
