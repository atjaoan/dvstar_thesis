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

def dvstar_build(threshold: float, min_count: int, max_depth: int):
    args = (
        "find",
        "./data/sequences_split_files",
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
        "./data/VLMCs/{}.bintree", # Fix the output path 
        ";"
    )
    subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def calculate_distances(dist_func: str) -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances", 
        "-p",
        cwd / "data/small_test/",
        "-n",
        dist_func 
    )
    res = subprocess.run(args, capture_output=True, text=True)

    return res

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
    res = subprocess.run(args, capture_output=True, text=True)

    return res

def save_to_csv(res: subprocess.CompletedProcess, csv_path: Path, dist_func: str):
    new_line_separated_attr = res.stderr.split('\n')[3:-8]

    print(res.stderr)

    data = [get_git_commit_version(), dist_func]
    columns = ["Repo Version", "distance_function"]
    
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

        data.extend([right_value, count])
        columns.extend([attribute, attribute + "_count"])

    new_line_separated_timings = res.stderr.split('\n')[-8:-3]
    for line in new_line_separated_timings:
        if len(line) < 1:
            continue
        split_line = line.lstrip().split(' ')
        if split_line[-1] == "elapsed":
            split_line[-1] = "elapsed time"
        data.append(split_line[0])
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
    timing_results = calculate_distances(dist_func.value)

    save_to_csv(timing_results, cwd / "tmp/benchmarks/test.csv", dist_func)

@app.command()
def record():
    dvstar_cmp_mem()

@app.command()
def build(threshold: float = 3.9075, min_count: int = 10, max_depth: int = 9):
    dvstar_build(threshold, min_count, max_depth)

if __name__ == "__main__":
    app()