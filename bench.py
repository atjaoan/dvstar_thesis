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

class VLMC_Container(str, Enum):
    vlmc_vector = "vector"
    vlmc_indexing = "indexing"
    vlmc_sorted_vector = "sorted-vector"
    vlmc_b_tree = "b-tree"
    vlmc_hashmap = "hashmap"
    vlmc_combo = "combo"
    vlmc_veb = "veb"
    vlmc_ey = "ey"
    vlmc_b_tree_alt = "alt-btree"
    vlmc_sorted_search = "sorted-search"

def get_bintree_name(genome_path: str, threshold: float, min_count: int, max_depth: int):
    return os.path.splitext(genome_path)[0] + f"_{threshold}_{min_count}_{max_depth}.bintree"

def get_parameter_from_bintree(bintree: str) -> tuple[float, int, int]:
    str_split = bintree.split('_')
    str_split.reverse()
    try:  
        return float(str_split[2]), int(str_split[1]), int(str_split[0].split('.')[0]) 
    except: 
        return 3.9075, 10, 9

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

def our_calculate_distances(dist_func: str, set_size: int, genome_path_fst: str, genome_path_snd, vlmc_container: str, nr_cores: int, background_order: int, mode: str) -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses,L1-dcache-load-misses,L1-dcache-loads,L1-dcache-stores,l2_rqsts.miss,LLC-load-misses,LLC-loads,LLC-prefetch-misses,LLC-store-misses,LLC-stores,cycle_activity.stalls_l3_miss",
        cwd / "build/dist", 
        "-p",
        cwd / genome_path_fst,
        "-s",
        cwd / genome_path_snd,
        "--function",
        dist_func,
        "-v",
        vlmc_container,
        "-n",
        str(nr_cores),
        "-b",
        str(background_order),
        "-m",
        mode
    )

    return subprocess.run(args, capture_output=True, text=True)

def our_calculate_distances_minimal(dist_func: str, set_size: int, genome_path_fst: str, genome_path_snd, vlmc_container: str, nr_cores: int, background_order: int, mode: str) -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "build/dist", 
        "-p",
        cwd / genome_path_fst,
        "-s",
        cwd / genome_path_snd,
        "--function",
        dist_func,
        "-v",
        vlmc_container,
        "-n",
        str(nr_cores),
        "-b",
        str(background_order),
        "-m",
        mode
    )

    return subprocess.run(args, capture_output=True, text=True)

def run_dev(nr_elements: int) -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-r 10",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "build/time_benchmark",
        str(nr_elements)
    )

    return subprocess.run(args, capture_output=True, text=True)

def cache_ref_to_vec_size(size: int) -> subprocess.CompletedProcess:
    args = (
        "perf",
        "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "build/dev", 
        str(size)
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
        "dvstar"
    )
    return subprocess.run(args, capture_output=True, text=True)

def save_to_csv(res: subprocess.CompletedProcess, csv_path: Path, dist_func: str, 
                set_size: int, threshold: float, min_count: int, max_depth: int, implementation: str):
    new_line_separated_attr = res.stderr.split('\n')[3:-8]

    print("For implementation -> " + implementation + "\n")
    print(res.stderr)
    return

    data = [get_git_commit_version(), implementation, dist_func, set_size, threshold, min_count, max_depth]
    columns = ["repo_version", "implementation", "distance_function", "set_size", "threshold", "min_count", "max_depth"]
    
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

def save_to_csv_dev(res: subprocess.CompletedProcess, csv_path: Path, nr_elements: int):
    new_line_separated_attr = res.stderr.split('\n')

    data = [get_git_commit_version(), nr_elements, res.stdout.split('\n')[0]]
    columns = ["repo_version", "nr_elements", "size_of_vector"]
    
    for line in new_line_separated_attr:
        split_line = line.split('#')
        if len(split_line) < 2:
            continue 

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

    new_line_separated_timings = res.stderr.split('\n')
    for line in new_line_separated_timings:
        if len(line) < 1:
            continue
        split_line = line.lstrip().split(' ')
        # print(split_line)        
        if split_line[5] == "elapsed":
            data.append(float(split_line[0].replace(",", ".")))
            columns.append("elapsed_time")
            break

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
    try: commit = res.stdout.split('\n')[0].split(' ')[1][0:7]
    except:
        try: commit = get_commit_from_file()
        except: 
            print("Failed to get commit hash, using empty string")
            return "" 
    return commit 

def get_commit_from_file():
    with open('current_commit.txt') as f:
        hash = f.readline()
        return hash[0:7]

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
        "dvstar"
    )
    subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    args = (
        "perf script | ./submodules/FlameGraph/stackcollapse-perf.pl | ./submodules/FlameGraph/flamegraph.pl > test.svg"
        )
    subprocess.run(args, shell=True)

    return

@app.command()
def dev_stat():
    for i in range(4, 22):
        print(f"Doing with {2**i}")
        res = run_dev(2**i)
        save_to_csv_dev(res, cwd / "csv_results/dev_test_6.csv", i)

@app.command()
def stat(set_size: int = -1, dist_func: Distance_Function = Distance_Function.dvstar, 
        genome_path_fst: str = "data/human_VLMCs", genome_path_snd: str = "data/human_VLMCs", background_order: int = 0):
    timing_results = calculate_distances(dist_func.value, -1, genome_path_fst, genome_path_snd, background_order)

    th, min, max = get_parameter_from_bintree(os.listdir(genome_path_fst)[0])

    save_to_csv(timing_results, cwd / "tmp/benchmarks/test.csv", dist_func, set_size, th, min, max, "Old")

@app.command()
def stat_new(set_size: int = -1, dist_func: Distance_Function = Distance_Function.dvstar, 
        genome_path_fst: str = "data/human_VLMCs", genome_path_snd: str = "data/human_VLMCs", vlmc_container: VLMC_Container = VLMC_Container.vlmc_combo, nr_cores: int = 1,
        background_order: int = 0, mode: str = "compare"):
    timing_results = our_calculate_distances(dist_func.value, set_size, genome_path_fst, genome_path_snd, vlmc_container.value, nr_cores, background_order, mode)

    th, min, max = get_parameter_from_bintree(os.listdir(genome_path_fst)[0])

    save_to_csv(timing_results, cwd / "tmp/benchmarks/test.csv", dist_func, set_size, th, min, max, "Opt-" + vlmc_container.value)

@app.command()
def stat_minimal(set_size: int = -1, dist_func: Distance_Function = Distance_Function.dvstar, 
        genome_path_fst: str = "data/human_VLMCs", genome_path_snd: str = "data/human_VLMCs", vlmc_container: VLMC_Container = VLMC_Container.vlmc_combo, nr_cores: int = 1,
        background_order: int = 0, mode: str = "compare"):
    timing_results = our_calculate_distances_minimal(dist_func.value, set_size, genome_path_fst, genome_path_snd, vlmc_container.value, nr_cores, background_order, mode)

    th, min, max = get_parameter_from_bintree(os.listdir(genome_path_fst)[0])

    save_to_csv(timing_results, cwd / "tmp/benchmarks/test.csv", dist_func, set_size, th, min, max, "Opt-" + vlmc_container.value)


@app.command()
def benchmark():
    background_order = 0
    genome_path_fst = "data/benchmarking/turkey/mega"
    genome_path_snd = "data/benchmarking/corn/mega"
    # genome_path_fst = "data/benchmarking/human/large"
    # genome_path_snd = "data/benchmarking/human/large"
    # stat(-1, Distance_Function.dvstar, genome_path_fst, genome_path_snd, background_order)
    stat_new(-1, Distance_Function.dvstar, genome_path_fst, genome_path_snd, VLMC_Container.vlmc_sorted_search, 8, background_order)
    # stat_new(-1, Distance_Function.dvstar, genome_path_fst, genome_path_snd, VLMC_Container.vlmc_ey, 8, background_order)
    # stat_new(-1, Distance_Function.dvstar, genome_path_fst, genome_path_snd, VLMC_Container.vlmc_veb, 8, background_order)
    stat_new(-1, Distance_Function.dvstar, genome_path_fst, genome_path_snd, VLMC_Container.vlmc_sorted_search, 8, background_order, "kmer-major")
    ## stat_new(-1, Distance_Function.dvstar, genome_path, VLMC_Container.vlmc_combo, 8, background_order)
    # stat_new(-1, Distance_Function.dvstar, genome_path, VLMC_Container.vlmc_combo, 8, background_order, "kmer-major")
    # stat_new(-1, Distance_Function.dvstar, genome_path_fst, genome_path_snd, VLMC_Container.vlmc_hashmap, 8)

@app.command()
def benchmarkmin():
    genome_path_fst = "data/benchmarking/human"
    print("#######################################")
    print("##                                   ##")
    print("##               HUMAN               ##")
    print("##                                   ##")
    print("#######################################")
    for folders in os.listdir(genome_path_fst):
        main_path = os.path.join(cwd, genome_path_fst)
        path = os.path.join(main_path, folders)
        stat_folder(path)
    genome_path_fst = "data/benchmarking/ecoli"
    print("#######################################")
    print("##                                   ##")
    print("##               ECOLI               ##")
    print("##                                   ##")
    print("#######################################")
    for folders in os.listdir(genome_path_fst):
        main_path = os.path.join(cwd, genome_path_fst)
        path = os.path.join(main_path, folders)
        stat_folder(path)

def stat_folder(path: str):
    background_order = 0
    stat(-1, Distance_Function.dvstar, path, path, background_order)
    stat_minimal(-1, Distance_Function.dvstar, path, path, VLMC_Container.vlmc_sorted_vector, 8, background_order)
    stat_minimal(-1, Distance_Function.dvstar, path, path, VLMC_Container.vlmc_ey, 8, background_order)
    stat_minimal(-1, Distance_Function.dvstar, path, path, VLMC_Container.vlmc_veb, 8, background_order)
    stat_minimal(-1, Distance_Function.dvstar, path, path, VLMC_Container.vlmc_b_tree, 8, background_order)
    stat_minimal(-1, Distance_Function.dvstar, path, path, VLMC_Container.vlmc_b_tree_alt, 8, background_order)
    stat_minimal(-1, Distance_Function.dvstar, path, path, VLMC_Container.vlmc_hashmap, 8)

@app.command()
def vecsize():
    r1 = cache_ref_to_vec_size(500)
    r2 = cache_ref_to_vec_size(5000)
    r3 = cache_ref_to_vec_size(50000)
    r4 = cache_ref_to_vec_size(500000)
    #r5 = cache_ref_to_vec_size(5000000)
    #r6 = cache_ref_to_vec_size(50000000)
    #r7 = cache_ref_to_vec_size(1000000000)
    print(r1.stderr)
    print(r1.stdout)
    print("-------------------------")
    print(r2.stderr)
    print(r2.stdout)
    print("-------------------------")
    print(r3.stderr)
    print(r3.stdout)
    print("-------------------------")
    print(r4.stderr)
    print(r4.stdout)
    #print("-------------------------")
    #print(r5.stderr)
    #print(r5.stdout)
    #print("-------------------------")
    #print(r6.stderr)
    #print(r6.stdout)
    #print("-------------------------")
    #print(r7.stderr)
    #print(r7.stdout)

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
def cache(dist_func: Distance_Function = Distance_Function.dvstar,
          genome_path: str = "data/human_VLMCs"):
    for size in range(50, 14500, 1450):
        stat(set_size=size, dist_func=dist_func, genome_path=genome_path)

if __name__ == "__main__":
    app()