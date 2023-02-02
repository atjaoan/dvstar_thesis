import typer
import subprocess
from pathlib import Path
from typing import Final
import re
import pandas as pd 
import os 
from enum import Enum
import numpy as np

app = typer.Typer()

K: Final[int] = 9

cwd = Path(__file__).parent

class Distance_Function(str, Enum):
    d2 = "d2"
    d2star = "d2star"
    dvstar = "dvstar"
    nearest_dvstar = "nearest-dvstar"
    penalized_dvstar = "penalized-dvstar"
    kl = "kl"
    kl_both = "kl-both"
    nll = "nll"
    nll_background = "nll-background"
    cv = "cv"
    cv_estimation = "cv-estimation"

class Build_Parameter(str, Enum):
    threshold = "threshold"
    min_count = "min-count"
    max_depth = "max-depth"

def get_bintree_name(genome_path: str, threshold: float, min_count: int, max_depth: int):
    return os.path.splitext(genome_path)[0] + f"_{threshold}_{min_count}_{max_depth}.bintree"

def get_parameter_from_bintree(bintree: str) -> tuple[float, int, int]:
    str_split = bintree.split('_')
    str_split.reverse()
    return float(str_split[2]), int(str_split[1]), int(str_split[0].split('.')[0]) 

def dvstar_build(genome_path: Path, out_path: Path, threshold: float, min_count: int, max_depth: int):
    for genome in os.listdir(genome_path):
        args = (
            "./build/dvstar",
            "--mode",
            "build",
            "--threshold",
            str(threshold),
            "--min-count",
            str(min_count),
            "--max-depth",
            str(max_depth),
            "--fasta-path",
            genome_path / genome,
            "--out-path",
            out_path / get_bintree_name(genome, threshold, min_count, max_depth)
        )
        subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def calculate_distances(dist_func: str, set_size: int, genome_path: str) -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances", 
        "-p",
        cwd / genome_path,
        "-n",
        dist_func,
        "-a",
        str(set_size)
    )

    return subprocess.run(args, capture_output=True, text=True)

def calculate_distances_only_oh() -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances-oh", 
        "-p",
        cwd / "data/small_test/",
        "-n",
        "d2"
    )
    return subprocess.run(args, capture_output=True, text=True)

def save_to_csv(res: subprocess.CompletedProcess, csv_path: Path, dist_func: str, 
                set_size: int, threshold: float, min_count: int, max_depth: int):
    new_line_separated_attr = res.stderr.split('\n')[3:-8]

    print(res.stderr)

    data = [get_git_commit_version(), dist_func, set_size, threshold, min_count, max_depth]
    columns = ["repo_version", "distance_function", "set_size", "threshold", "min_count", "max_depth"]
    
    for line in new_line_separated_attr:
        split_line = line.split('#')

        right_value = split_line[1].lstrip().split(' ')[0]
        left_line = split_line[0].strip()

        if "msec" in left_line:
            count_and_attribute = left_line.split('msec')
        else: 
            count_and_attribute = re.split(r"\s{2,}", left_line)

        attribute  = count_and_attribute[1].replace(':u', '').strip().replace('-', '_')
        count = count_and_attribute[0]
        
        # Remove % if it exists and replace commas with dot
        right_value = right_value.strip('%').replace(',', '.')
        
        # For counts to be made into ints, skip space separator
        count = count.replace('\u202f', '').replace(',', '.').rstrip()

        data.extend([float(right_value), int(float(count))])
        columns.extend([attribute, attribute + "_count"])

    new_line_separated_timings = res.stderr.split('\n')[-8:-3]
    for line in new_line_separated_timings:
        if len(line) < 1:
            continue
        split_line = line.lstrip().split(' ')
        if split_line[-1] == "elapsed":
            split_line[-1] = "elapsed_time"
        data.append(float(split_line[0].replace(",", ".")))
        columns.append(split_line[-1])

    if not os.path.exists(csv_path):
        df = pd.DataFrame(columns=columns)
    else:
        df = pd.read_csv(csv_path, dtype=str)

    df.loc[len(df)] = data
    df.to_csv(csv_path, index=False)

def get_git_commit_version():
    args = (
        "git",
        "log",
        "-1"
    )
    res = subprocess.run(args, capture_output=True, text=True)
    commit = res.stdout.split('\n')[0].split(' ')[1][0:7]
    return commit 

def dvstar_cmp_mem():
    args = (
        "perf",
        "record",
        "--call-graph",
        "dwarf",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances",
        "-p",
        cwd / "data/small_test/",
        "-n",
        "d2"
    )
    subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    args = (
        "perf script | ./submodules/FlameGraph/stackcollapse-perf.pl | ./submodules/FlameGraph/flamegraph.pl > test.svg"
        )
    subprocess.run(args, shell=True)

    return

@app.command()
def stat(set_size: int = -1, dist_func: Distance_Function = Distance_Function.dvstar, 
        genome_path: str = "data/human_VLMCs"):
    timing_results = calculate_distances(dist_func.value, -1, genome_path)

    th, min, max = get_parameter_from_bintree(os.listdir(genome_path)[0])

    save_to_csv(timing_results, cwd / "tmp/benchmarks/test.csv", dist_func, set_size, th, min, max)

@app.command()
def record():
    dvstar_cmp_mem()

@app.command()
def build(threshold: float = 3.9075, min_count: int = 10, max_depth: int = 9, 
        genome_path: str ="./data/test", out_path: str ="./data/test_VLMCs"):
    genome_path = Path(genome_path)
    out_path    = Path(out_path)
    dvstar_build(genome_path, out_path, threshold, min_count, max_depth)

@app.command()
def human_build(threshold: float = 3.9075, min_count: int = 10, max_depth: int = 9):
    build(threshold, min_count, max_depth, "./data/human_genome_split_files", "./data/human_VLMCs")

@app.command()
def cache(dist_func: Distance_Function = Distance_Function.dvstar,
          genome_path: str = "data/human_VLMCs"):
    for size in range(50, 14500, 1450):
        stat(set_size=size, dist_func=dist_func, genome_path=genome_path)

@app.command()
def vlmc_size_benchmark(parameter: Build_Parameter, genome_path: str = "./data/test"):
    vals = []
    if parameter == Build_Parameter.threshold:
        for x in np.arange(2.0, 4.0, 0.25):
            out_path=f"./data/vlmc_size_benchmark/{parameter}_{x}"
            build(threshold=x, genome_path=genome_path, out_path=out_path)
            vals.append(x)
    elif parameter == Build_Parameter.min_count:
        for x in range(4, 18, 2):
            out_path=f"./data/vlmc_size_benchmark/{parameter}_{x}"
            build(min_count=x, genome_path=genome_path, out_path=out_path)
            vals.append(x)
    else: 
        for x in range(4, 18, 2):
            out_path=f"./data/vlmc_size_benchmark/{parameter}_{x}"
            build(max_depth=x, genome_path=genome_path, out_path=out_path)
            vals.append(x)     

    for x in vals:
        print(f"With {parameter} set to {x}")
        stat(genome_path=f"./data/vlmc_size_benchmark/{parameter}_{x}")


if __name__ == "__main__":
    app()