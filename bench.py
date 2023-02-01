import typer
import subprocess
from pathlib import Path
from typing import Final
import re
import pandas as pd 
import os 
from enum import Enum

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

def dvstar_build(genome_path: Path, threshold: float, min_count: int, max_depth: int):
    args = (
        "find",
        genome_path,
        "-name",
        "*.fasta",
        "-exec",
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
        "{}",
        "--out-path",
        "./data/human/{}.bintree", # Fix the output path 
        ";"
    )
    subprocess.run(args)#, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL

def calculate_distances(dist_func: str, set_size: int) -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances", 
        "-p",
        cwd / "data/VLMCs",
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

def save_to_csv(res: subprocess.CompletedProcess, csv_path: Path, dist_func: str, set_size: int):
    new_line_separated_attr = res.stderr.split('\n')[3:-8]

    print(res.stderr)

    data = [get_git_commit_version(), dist_func, set_size]
    columns = ["Repo Version", "distance_function", "Set size"]
    
    for line in new_line_separated_attr:
        split_line = line.split('#')

        right_value = split_line[1].lstrip().split(' ')[0]
        left_line = split_line[0].strip()

        if "msec" in left_line:
            count_and_attribute = left_line.split('msec')
        else: 
            count_and_attribute = re.split(r"\s{2,}", left_line)

        attribute  = count_and_attribute[1].replace(":u", "")
        count = count_and_attribute[0]
        
        # Remove % if it exists and replace commas with dot
        right_value = right_value.strip("%").replace(",", ".")
        
        # For counts to be made into ints, skip space separator
        count = count.replace("\u202f", "").replace(",", ".").rstrip()

        data.extend([float(right_value), int(float(count))])
        columns.extend([attribute, attribute + "_count"])

    new_line_separated_timings = res.stderr.split('\n')[-8:-3]
    for line in new_line_separated_timings:
        if len(line) < 1:
            continue
        split_line = line.lstrip().split(' ')
        if split_line[-1] == "elapsed":
            split_line[-1] = "elapsed time"
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
def stat(dist_func: Distance_Function = Distance_Function.dvstar):
    timing_results = calculate_distances(dist_func.value, -1)

    save_to_csv(timing_results, cwd / "tmp/benchmarks/test.csv", dist_func)

@app.command()
def record():
    dvstar_cmp_mem()

@app.command()
def build(threshold: float = 3.9075, min_count: int = 10, max_depth: int = 9):
    dvstar_build("./data/sequences_split_files", threshold, min_count, max_depth)

@app.command()
def human_build(threshold: float = 3.9075, min_count: int = 10, max_depth: int = 9):
    dvstar_build("./data/human_genome_split_files", threshold, min_count, max_depth)

@app.command()
def cache(dist_func: Distance_Function = Distance_Function.dvstar):
    for size in range(50, 14500, 1450):
        timing_results = calculate_distances(dist_func, size)

        save_to_csv(timing_results, cwd / "tmp/benchmarks/test.csv", dist_func, size)

if __name__ == "__main__":
    app()