import typer
import subprocess
from pathlib import Path
import os 
import numpy as np 

app = typer.Typer()

cwd = Path(__file__).parent

def get_bintree_name(genome_path: str, threshold: float, min_count: int, max_depth: int):
    return os.path.splitext(genome_path)[0] + f"_{threshold}_{min_count}_{max_depth}.bintree"

#####################################################
# Builds VLMCs from sequences.                      #
# Dvstars projects implementation is used.          #
# threshold = 3.9075, min_count = 10, max_depth = 9 # 
#####################################################
def dvstar_build(genome_path: Path, out_path: Path, threshold: float, min_count: int, max_depth: int):
    for genome in os.listdir(genome_path):
        args = (
            "./build/dvstar",
            "--mode", "build",
            "--threshold", str(threshold),
            "--min-count", str(min_count),
            "--max-depth", str(max_depth),
            "--fasta-path", genome_path / genome,
            "--out-path", out_path / get_bintree_name(genome, threshold, min_count, max_depth)
        )
        subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

#####################################
# Build data set for comparing file #
# size to input parameters.         #
#####################################
def build_parameter_test():
    genome_path = cwd / "data/test"
    out_path = cwd / "data/benchmarking/parameter_test"
    for threshold in np.arange(0, 4, 0.25):
        for min_count in range(0, 15, 1):
            for max_depth in np.arange(0, 15, 1):
                dvstar_build(genome_path, out_path, threshold, min_count, max_depth) 

#########################
# Build human data set. #
#########################
def build_human():
    genome_path = cwd / "data/human_genome_split_files"
    out_path = cwd / "data/benchmarking/human"
    print("Building small vlmcs...")
    dvstar_build(genome_path, out_path / "small", 3.9075, 9, 6)
    print("Building medium vlmcs...")
    dvstar_build(genome_path, out_path / "medium", 3, 6, 8)
    print("Building large vlmcs...")
    dvstar_build(genome_path, out_path / "large", 2, 3, 10)

##########################
# Build E-coli data set. #
##########################
def build_ecoli():
    genome_path = cwd / "data/sequences_split_files"
    out_path = cwd / "data/benchmarking/ecoli"
    print("Building small vlmcs...")
    dvstar_build(genome_path, out_path / "small", 3.9075, 9, 6)
    print("Building medium vlmcs...")
    dvstar_build(genome_path, out_path / "medium", 3, 6, 8)
    print("Building large vlmcs...")
    dvstar_build(genome_path, out_path / "large", 2, 3, 10)

@app.command()
def build():
    ## print("Building Human Sequences...")
    ## build_human()
    print("Building E-coli Sequences...")
    build_ecoli()


if __name__ == "__main__":
    app()