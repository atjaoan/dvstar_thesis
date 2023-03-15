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
    
def get_csv_name():
    now = datetime.now()
    return "csv_results/" + now.strftime("%m_%d_%H_%M") + ".csv"

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
def PstClassifierSeqan(set_size: int, genome_path: str, background_order: int) -> subprocess.CompletedProcess:
    args = (
        "perf", "stat",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances", 
        "-p", cwd / genome_path,
        "-n", "dvstar",
        "-a", str(set_size),
        "-b", str(background_order)
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

#####################################
# Function for running all tests on #
# the human dataset.                #
#####################################
def human_benchmarking():
    dir_path = cwd / "data/benchmarking/human"
    csv_filename = get_csv_name()
    nb_files = count_nb_files(dir_path / "small")

    th_small, min_small, max_small = get_parameter_from_bintree(os.listdir(dir_path / "small")[0])
    th_medium, min_medium, max_medium = get_parameter_from_bintree(os.listdir(dir_path / "medium")[0])
    th_large, min_large, max_large = get_parameter_from_bintree(os.listdir(dir_path / "large")[0])

    while(nb_files > 2):
        print("Benchmarking with " + str(nb_files) + " VLMCs...")
        print("Benchmarking small PstClassifierSeqan.")
        res_small  = PstClassifierSeqan(nb_files, dir_path / "small", 0)
        print("Benchmarking medium PstClassifierSeqan.")
        res_medium = PstClassifierSeqan(nb_files, dir_path / "medium", 0)
        print("Benchmarking large PstClassifierSeqan.")
        res_large  = PstClassifierSeqan(nb_files, dir_path / "large", 0)

        save_to_csv(res_small, cwd / csv_filename, "small", nb_files, th_small, min_small, max_small, "PstClassifierSeqan", 8)
        save_to_csv(res_medium, cwd / csv_filename, "medium", nb_files, th_medium, min_medium, max_medium, "PstClassifierSeqan", 8)
        save_to_csv(res_large, cwd / csv_filename, "large", nb_files, th_large, min_large, max_large, "PstClassifierSeqan", 8)

        for imp in ['sorted-vector']:
            print("Benchmarking small " + imp + ".")
            res_our_small  = calculate_distance(nb_files, dir_path / "small", imp, 8, 0)
            print("Benchmarking medium " + imp + ".")
            res_our_medium = calculate_distance(nb_files, dir_path / "medium", imp, 8, 0)
            print("Benchmarking large " + imp + ".")
            res_our_large  = calculate_distance(nb_files, dir_path / "large", imp, 8, 0)

            save_to_csv(res_our_small, cwd / csv_filename, "small", nb_files, th_small, min_small, max_small, imp, 8)
            save_to_csv(res_our_medium, cwd / csv_filename, "medium", nb_files, th_medium, min_medium, max_medium, imp, 8)
            save_to_csv(res_our_large, cwd / csv_filename, "large", nb_files, th_large, min_large, max_large, imp, 8)

        nb_files = int(nb_files / 2)


@app.command()
def benchmark():
    human_benchmarking()

if __name__ == "__main__":
    app()