import typer
import subprocess
from pathlib import Path
from typing import Final
import time
import re
import pandas as pd 
import os 

app = typer.Typer()

K: Final[int] = 9

cwd = Path(__file__).parent

def dvstar_build():
    args = (
        "find",
        "./data/sequences_split_files",
        "-name",
        "*.fasta",
        "-exec",
        "./build/dvstar",
        "--mode",
        "build",
        "--min-count",
        "10",
        "--max-depth",
        "9",
        "--fasta-path",
        "{}",
        "--out-path",
        "./data/VLMCs/{}.bintree", # Fix the output path 
        ";"
    )
    subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def dvstar_cmp(path_1: Path, path_2: Path) -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "build/dvstar",
        "--mode",
        "5",
        "--in-path",
        path_1,
        "--to-path",
        path_2
    )
    res = subprocess.run(args, capture_output=True, text=True)

    return res

def calculate_distances() -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances", # Hard coded path to calculate-distances might need to be changed
        "-p",
        cwd / "data/small_test/",
        "-n",
        "d2"
    )
    res = subprocess.run(args, capture_output=True, text=True)

    return res

def calculate_distances_only_oh() -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances-oh", # Hard coded path to calculate-distances might need to be changed
        "-p",
        cwd / "data/small_test/",
        "-n",
        "d2"
    )
    res = subprocess.run(args, capture_output=True, text=True)

    return res

def save_to_csv(res: subprocess.CompletedProcess, csv_path: Path):
    new_line_separated_attr = res.stderr.split('\n')[3:-8]

    print(res.stderr)

    data = [get_git_commit_version()]
    columns = ["Repo Version"]
    
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

def dvstar_cmp_mem(path_1: Path, path_2: Path):
    args = (
        "perf",
        "record",
        "-F 99",
        cwd / "build/dvstar",
        "--mode",
        "5",
        "--in-path",
        path_1,
        "--to-path",
        path_2
    )
    subprocess.run(args)

    args = (
        "perf script |",
        cwd / "submodules" / "FlameGraph" / "./stackcollapse-perf.pl |",
        cwd / "submodules" / "FlameGraph" / "./flamegraph.pl > call_stack.svg"
    )
    return

@app.command()
def stat():
    timing_results = calculate_distances_only_oh()

    save_to_csv(timing_results, cwd / "tmp/benchmarks/test.csv")

@app.command()
def build():
    dvstar_build()

if __name__ == "__main__":
    app()