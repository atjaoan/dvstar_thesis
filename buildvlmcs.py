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
            # "--in-or-out-of-core",
            # "external"
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

###########################
# Check if dir exists and # 
# pass promt if true.     #
###########################
def promt_dir_build(path : Path, size : str, name : str):
    if (os.path.isdir(path / size)):
        x = input("The " + size + " dataset for " + name + " already exists. Do you want to build it anyways? (Y/n) \n")
        while (x != "n") & (x != "Y"):
            x = input("Could not interpret input. Please write 'Y' or 'n' \n")
        if x == "n":
            return False
    print("Building " + size + " vlmcs...")
    return True 

#########################
# Build data set. #
#########################
def build_dataset(path_to_fasta : str, name : str):
    genome_path = cwd / path_to_fasta
    out_path = cwd / "data/benchmarking" / name
    if (os.path.isdir(genome_path)):
        if promt_dir_build(out_path, "small", name):
            dvstar_build(genome_path, out_path / "small", 3.9075, 9, 6)
        if promt_dir_build(out_path, "medium", name):
            dvstar_build(genome_path, out_path / "medium", 3.9075, 12, 10) # 3, 6, 8 old
        if promt_dir_build(out_path, "large", name):
            dvstar_build(genome_path, out_path / "large", 3.9075, 15, 15) # 2, 3, 10 old 
    else: 
        print("Could not find the directory for creating '" + name + "' benchmarking files.")
        print("The given directory was : '" + str(genome_path) + "'")

@app.command()
def build():
    print("Building Human Sequences...")
    build_dataset("data/human_genome_split_files", "human")
    print("Building E-coli Sequences...")
    build_dataset("data/ecoli_split_files", "ecoli")
    print("Building Turkey Sequences...")
    build_dataset("data/turkey_fasta_files", "turkey")
    print("Building Corn Sequences...")
    build_dataset("data/corn_fasta_files", "corn")


if __name__ == "__main__":
    app()