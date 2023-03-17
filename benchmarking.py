import typer
import subprocess
from pathlib import Path
from typing import Final
import re
import pandas as pd 
import os 
from enum import Enum
import numpy as np
from datetime import datetime

app = typer.Typer()

cwd = Path(__file__).parent

#############################
# VLMC container available  #
# to use for benchmarking.  #
#############################
class VLMC_Container(str, Enum):
    vlmc_vector = "vector"
    vlmc_indexing = "indexing"
    vlmc_sorted_vector = "sorted-vector"
    vlmc_b_tree = "b-tree"
    vlmc_hashmap = "hashmap"
    vlmc_combo = "combo"

##############################
# Get current git version to #
# save in csv file.          #
##############################
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
    
def get_csv_name(dataset: str, extra: str = ""):
    dset = dataset.split("/")
    now = datetime.now()
    return "csv_results/" + extra + dset[len(dset) - 1] + "_" + now.strftime("%m_%d_%H_%M") + ".csv"

def count_nb_files(dir: str):
    count = 0
    for file in os.listdir(dir):
        count += 1
    return count 

##############################
# Get build parameters from  #
# filename.                  #
##############################
def get_parameter_from_bintree(bintree: str) -> tuple[float, int, int]:
    str_split = bintree.split('_')
    str_split.reverse()
    try:  
        return float(str_split[2]), int(str_split[1]), int(str_split[0].split('.')[0]) 
    except: 
        return 3.9075, 10, 9

####################################
# Call to our implementation of    #
# distance calculation between     #
# one or two directories of VLMCs. #
####################################     
def calculate_distance(set_size: int, genome_path: str, vlmc_container: str, nr_cores: int, background_order: int) -> subprocess.CompletedProcess:
    args = (
        "perf", "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "build/dist", 
        "-p", cwd / genome_path,
        "-s", cwd / genome_path,
        "-v", vlmc_container,
        "-n", str(nr_cores),
        "-b", str(background_order),
        "-a", str(set_size)
    )
    return subprocess.run(args, capture_output=True, text=True)

######################################
# Call to the old implementation of  #
# distance calculation between VLMCs #
# from PstClassifierSeqan submodule. #
######################################
def PstClassifierSeqan(set_size: int, genome_path: str, out_path: str, background_order: int) -> subprocess.CompletedProcess:
    args = (
        "perf", "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances", 
        "-p", cwd / genome_path,
        "-n", "dvstar",
        "-a", str(set_size),
        "-b", str(background_order),
        "-s", cwd / out_path
    )
    return subprocess.run(args, capture_output=True, text=True)

######################################
# Save output from perf to csv file. #
######################################
def save_to_csv(res: subprocess.CompletedProcess, csv_path: Path, vlmc_size: str, 
                set_size: int, threshold: float, min_count: int, max_depth: int, implementation: str, nr_cores_used: int):
    new_line_separated_attr = res.stderr.split('\n')[3:-8]

    print("For implementation -> " + implementation + "\n")
    print(res.stderr)

    data = [get_git_commit_version(), implementation, vlmc_size, set_size, threshold, min_count, max_depth, nr_cores_used]
    columns = ["repo_version", "implementation", "vlmc_size", "set_size", "threshold", "min_count", "max_depth", "nr_cores_used"]
    
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

######################################
# Function for benchmarking time and # 
# cache-misses on dataset.           #
######################################
def normal_benchmaking(dataset: str, implementation: str):
    dir_path = cwd / dataset
    csv_filename = get_csv_name(dataset)
    print(csv_filename)
    nb_files = count_nb_files(dir_path / "small")

    nb_files = int(nb_files / 2)

    th_small, min_small, max_small = get_parameter_from_bintree(os.listdir(dir_path / "small")[0])
    th_medium, min_medium, max_medium = get_parameter_from_bintree(os.listdir(dir_path / "medium")[0])
    th_large, min_large, max_large = get_parameter_from_bintree(os.listdir(dir_path / "large")[0])

    while(nb_files > 2):
        print("Benchmarking with " + str(nb_files) + " VLMCs...")
        print("Benchmarking small " + implementation + ".")
        res_our_small  = calculate_distance(nb_files, dir_path / "small", implementation, 8, 0)
        print("Benchmarking medium " + implementation + ".")
        res_our_medium = calculate_distance(nb_files, dir_path / "medium", implementation, 8, 0)
        print("Benchmarking large " + implementation + ".")
        res_our_large  = calculate_distance(nb_files, dir_path / "large", implementation, 8, 0)

        save_to_csv(res_our_small, cwd / csv_filename, "small", nb_files, th_small, min_small, max_small, implementation, 8)
        save_to_csv(res_our_medium, cwd / csv_filename, "medium", nb_files, th_medium, min_medium, max_medium, implementation, 8)
        save_to_csv(res_our_large, cwd / csv_filename, "large", nb_files, th_large, min_large, max_large, implementation, 8)

        nb_files = int(nb_files / 2)

def Pst_normal_benchmaking(dataset: str):
    dir_path = cwd / dataset
    csv_filename = get_csv_name(dataset)
    print(csv_filename)
    nb_files = count_nb_files(dir_path / "small")

    nb_files = int(nb_files / 2)

    th_small, min_small, max_small = get_parameter_from_bintree(os.listdir(dir_path / "small")[0])
    th_medium, min_medium, max_medium = get_parameter_from_bintree(os.listdir(dir_path / "medium")[0])
    th_large, min_large, max_large = get_parameter_from_bintree(os.listdir(dir_path / "large")[0])

    while(nb_files > 2):
        print("Benchmarking small PstClassifierSeqan.")
        res_small  = PstClassifierSeqan(nb_files, dir_path / "small", 0)
        print("Benchmarking medium PstClassifierSeqan.")
        res_medium = PstClassifierSeqan(nb_files, dir_path / "medium", 0)
        print("Benchmarking large PstClassifierSeqan.")
        res_large  = PstClassifierSeqan(nb_files, dir_path / "large", 0)

        save_to_csv(res_small, cwd / csv_filename, "small", nb_files, th_small, min_small, max_small, "PstClassifierSeqan", 8)
        save_to_csv(res_medium, cwd / csv_filename, "medium", nb_files, th_medium, min_medium, max_medium, "PstClassifierSeqan", 8)
        save_to_csv(res_large, cwd / csv_filename, "large", nb_files, th_large, min_large, max_large, "PstClassifierSeqan", 8)
        nb_files = int(nb_files / 2)

#####################################
# Function for benchmarking degree  #
# of parallelization on dataset.    #  
#####################################
def parallelization_benchmark(dataset: str):
    dir_path = cwd / dataset
    csv_filename = get_csv_name(dataset, "parallelization_")
    print(csv_filename)
    nb_files = count_nb_files(dir_path / "small")

    nb_files = int(nb_files / 2)

    th_small, min_small, max_small = get_parameter_from_bintree(os.listdir(dir_path / "small")[0])
    th_medium, min_medium, max_medium = get_parameter_from_bintree(os.listdir(dir_path / "medium")[0])
    th_large, min_large, max_large = get_parameter_from_bintree(os.listdir(dir_path / "large")[0])
    for nr_cores in range(1, 9, 1):
        print("Benchmarking small with " + str(nr_cores) + " cores.")
        res_our_small  = calculate_distance(nb_files, dir_path / "small", 'sorted-vector', nr_cores, 0)
        print("Benchmarking medium with " + str(nr_cores) + " cores.")
        res_our_medium = calculate_distance(nb_files, dir_path / "medium", 'sorted-vector', nr_cores, 0)
        print("Benchmarking large with " + str(nr_cores) + " cores.")
        res_our_large  = calculate_distance(nb_files, dir_path / "large", 'sorted-vector', nr_cores, 0)

        save_to_csv(res_our_small, cwd / csv_filename, "small", nb_files, th_small, min_small, max_small, 'sorted-vector', nr_cores)
        save_to_csv(res_our_medium, cwd / csv_filename, "medium", nb_files, th_medium, min_medium, max_medium, 'sorted-vector', nr_cores)
        save_to_csv(res_our_large, cwd / csv_filename, "large", nb_files, th_large, min_large, max_large, 'sorted-vector', nr_cores)


import h5py

def compare_hdf5_files(vlmc_size: str, PstGroup: str, VLMCGroup: str):
    vlmc_file = h5py.File(cwd / ("hdf5_results/distances" + vlmc_size + ".hdf5"), 'r')
    pst_file  = h5py.File(cwd / ("hdf5_results/pst_distances" + vlmc_size + ".hdf5"), 'r')

    vlmc_group = vlmc_file.get('distances').get(VLMCGroup)
    pst_group  = pst_file.get('distances').get(PstGroup)

    print(list(pst_file.get('distances').keys()))

    error_tolerance = 1e-8

    for x in range(0, vlmc_group[0].size):
        for y in range(0, vlmc_group[0].size):
            if abs(vlmc_group[x][y] - pst_group[x][y]) > error_tolerance:
                print("Our = " + str(vlmc_group[x][y]) + " pst = " + str(pst_group[x][y]))
                print("Difference = " + str(vlmc_group[x][y] - pst_group[x][y]))

@app.command()
def benchmark():
    PstClassifierSeqan(-1, "data/test_VLMCs", "hdf5_results/pst_distances.hdf5", 0)
    compare_hdf5_files("", "dvstar-0", "sorted-vector-test_VLMCs-test_VLMCs")
    ## Pst_normal_benchmaking("data/benchmarking/ecoli")
    ## parallelization_benchmark("data/benchmarking/ecoli")

if __name__ == "__main__":
    app()