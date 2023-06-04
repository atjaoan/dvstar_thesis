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

class Build_Parameter(str, Enum):
    threshold = "threshold"
    min_count = "min-count"
    max_depth = "max-depth"

class VLMC_Container(str, Enum):
    vlmc_sorted_vector = "sorted-vector"
    vlmc_hashmap = "hashmap"
    vlmc_veb = "veb"
    vlmc_ey = "eytzinger"
    vlmc_b_tree_alt = "b-tree"
    vlmc_sorted_search = "sbs"

def calculate_distances(dist_func: str, set_size: int, genome_path_fst: str, genome_path_snd: str, background_order: int) -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances", 
        "-p",
        cwd / genome_path_fst,
        "-s",
        cwd / genome_path_snd,
        "-n",
        dist_func,
        "-a",
        str(set_size),
        "-b",
        str(background_order)
    )

    return subprocess.run(args, capture_output=True, text=True)

def our_calculate_distances(set_size: int, genome_path_fst: str, genome_path_snd, vlmc_container: str, nr_cores: int, background_order: int, mode: str) -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses,L1-dcache-load-misses,L1-dcache-loads,L1-dcache-stores,l2_rqsts.miss,LLC-load-misses,LLC-loads,LLC-prefetch-misses,LLC-store-misses,LLC-stores,cycle_activity.stalls_l3_miss",
        cwd / "build/dist", 
        "-p",
        cwd / genome_path_fst,
        "-s",
        cwd / genome_path_snd,
        "-v",
        vlmc_container,
        "-n",
        str(nr_cores),
        "-b",
        str(background_order),
        "-a",
        str(set_size)
    )

    return subprocess.run(args, capture_output=True, text=True)

@app.command()
def stat_new(set_size: int = -1, 
        genome_path_fst: str = "data/human_VLMCs", genome_path_snd: str = "data/human_VLMCs", vlmc_container: VLMC_Container = VLMC_Container.vlmc_sorted_search, nr_cores: int = 1,
        background_order: int = 0, mode: str = "compare"):
    timing_results = our_calculate_distances(set_size, genome_path_fst, genome_path_snd, vlmc_container.value, nr_cores, background_order, mode)
    print(timing_results.stderr)

@app.command()
def benchmark():
    background_order = 0
    genome_path_fst = "data/benchmarking/human/large"
    genome_path_snd = "data/benchmarking/human/large"
    # genome_path_fst = "data/benchmarking/human/large"
    # genome_path_snd = "data/benchmarking/human/large"
    # stat(-1, Distance_Function.dvstar, genome_path_fst, genome_path_snd, background_order)
    stat_new(-1, genome_path_fst, genome_path_snd, VLMC_Container.vlmc_sorted_vector, 8, background_order)
    stat_new(-1, genome_path_fst, genome_path_snd, VLMC_Container.vlmc_sorted_search, 8, background_order)
    stat_new(-1, genome_path_fst, genome_path_snd, VLMC_Container.vlmc_hashmap, 8, background_order)

if __name__ == "__main__":
    app()