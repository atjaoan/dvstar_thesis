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
import h5py

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
def calculate_distance(set_size: int, genome_path: str, out_path: str, vlmc_container: str, nr_cores: int, background_order: int) -> subprocess.CompletedProcess:
    args = (
        "perf", "stat",
        "-r",
        "5",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "build/dist", 
        "-p", cwd / genome_path,
        "-s", cwd / genome_path,
        "-v", vlmc_container,
        "-n", str(nr_cores),
        "-b", str(background_order),
        "-a", str(set_size),
        "-o", str(out_path)
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
        "-r",
        "5",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances", 
        "-p", cwd / genome_path,
        "-y", cwd / genome_path,
        "-n", "dvstar",
        "-a", str(set_size),
        "-b", str(background_order),
        "-s", cwd / out_path
    )
    return subprocess.run(args, capture_output=True, text=True)

######################################
# Save output from perf to csv file. #
######################################
def catch_and_save(res: subprocess.CompletedProcess, csv_path: Path, vlmc_size: str, 
                set_size: int, threshold: float, min_count: int, max_depth: int, 
                implementation: str, nr_cores_used: int):
    try: 
        save_to_csv(res, csv_path, vlmc_size, set_size, threshold, min_count, max_depth, implementation, nr_cores_used)
    except:
        print(res.stderr)
        raise Exception("Failed to write to csv")

def save_to_csv(res: subprocess.CompletedProcess, csv_path: Path, vlmc_size: str, 
                set_size: int, threshold: float, min_count: int, max_depth: int, 
                implementation: str, nr_cores_used: int, combo_init_size: int = 128):
    new_line_separated_attr = res.stderr.split('\n')[3:-4]

    print(f"Save implementation -> {implementation} to csv...")

    data = [get_git_commit_version(), implementation, vlmc_size, set_size, threshold, min_count, max_depth, nr_cores_used, combo_init_size]
    columns = ["repo_version", "implementation", "vlmc_size", "set_size", "threshold", "min_count", "max_depth", "nr_cores_used", "combo_init_size"]
    
    for line in new_line_separated_attr:
        split_line = line.split('#')

        divergence = split_line[1].split('(')[1].lstrip().split(' ')[1]
        if (split_line[1].split('(')[1].lstrip().split(' ')[1] == ''):
            divergence = split_line[1].split('(')[1].lstrip().split(' ')[2]

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

        # Remove % if it exists and replace commas with dot
        divergence = divergence.strip('%').replace(',', '.')
        
        # For counts to be made into ints, skip space separator
        count = count.replace('\u202f', '').replace(',', '.').rstrip()

        data.extend([float(right_value), int(float(count)), float(divergence)])
        columns.extend([attribute, attribute + "_count", attribute + "_divergence"])

    new_line_separated_timings = res.stderr.split('\n')[-3]
    split_line = new_line_separated_timings.split('+-')
    divergence = split_line[1].lstrip().split(' ')[0]
    data.extend([float(split_line[0].replace(",", ".")), float(divergence.replace(",", "."))])
    columns.extend(["elapsed_time", "elapsed_time_divergence"])

    if not os.path.exists(csv_path):
        df = pd.DataFrame(columns=columns)
    else:
        df = pd.read_csv(csv_path, dtype=str)

    df.loc[len(df)] = data
    df.to_csv(csv_path, index=False)

######################################
# Function for benchmarking time and # 
# cache-misses on dataset with       #
# Pst... implementation.             #
######################################
def Pst_normal_benchmaking(dataset: str):
    csv_filename = get_csv_name(dataset)
    print(csv_filename)
    nb_files = count_nb_files(cwd / dataset / "small")

    files_run = 2   

    th_small, min_small, max_small = get_parameter_from_bintree(os.listdir(cwd / dataset / "small")[0])
    th_medium, min_medium, max_medium = get_parameter_from_bintree(os.listdir(cwd / dataset / "medium")[0])
    th_large, min_large, max_large = get_parameter_from_bintree(os.listdir(cwd / dataset / "large")[0])

    out_path = "hdf5_results/PstClassifierSeqan/" + dataset + "/"
    
    os.makedirs(out_path, exist_ok=True)

    while((files_run < 10000) & (files_run < nb_files)):
        print("Benchmarking with " + str(files_run) + " VLMCs...")
        print("Small PstClassifierSeqan.")
        res_small  = PstClassifierSeqan(files_run, dataset + "/small", out_path + "small_" + str(files_run) + ".hdf5", 0)
        print("Medium PstClassifierSeqan.")
        res_medium = PstClassifierSeqan(files_run, dataset + "/medium", out_path + "medium_" + str(files_run) + ".hdf5", 0)
        print("Large PstClassifierSeqan.")
        res_large  = PstClassifierSeqan(files_run, dataset + "/large", out_path + "large_" + str(files_run) + ".hdf5", 0)

        catch_and_save(res_small, cwd / csv_filename, "small", files_run, th_small, min_small, max_small, "PstClassifierSeqan", 8)
        catch_and_save(res_medium, cwd / csv_filename, "medium", files_run, th_medium, min_medium, max_medium, "PstClassifierSeqan", 8)
        catch_and_save(res_large, cwd / csv_filename, "large", files_run, th_large, min_large, max_large, "PstClassifierSeqan", 8)
        files_run = files_run * 2

######################################
# Function for benchmarking time and # 
# cache-misses on dataset with new   #
# implementation.                    #
######################################
def normal_benchmaking(dataset: str, implementation: str):
    csv_filename = get_csv_name(dataset, implementation + "_")
    print(csv_filename)
    nb_files = count_nb_files(cwd / dataset / "small")

    th_small, min_small, max_small = get_parameter_from_bintree(os.listdir(cwd / dataset / "small")[0])
    th_medium, min_medium, max_medium = get_parameter_from_bintree(os.listdir(cwd / dataset / "medium")[0])
    th_large, min_large, max_large = get_parameter_from_bintree(os.listdir(cwd / dataset / "large")[0])

    files_run = 2

    out_path = "hdf5_results/Dvstar/"

    os.makedirs(out_path, exist_ok=True)

    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_small.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_small.hdf5")
    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")
    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_large.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_large.hdf5")

    while((files_run < 10000) & (files_run < nb_files)):
        print("Benchmarking with " + str(files_run) + " VLMCs...")
        print("Small " + implementation + ".")
        res_our_small  = calculate_distance(files_run, dataset + "/small", "hdf5_results/Dvstar/distances_small.hdf5", implementation, 8, 0)
        print("Medium " + implementation + ".")
        res_our_medium = calculate_distance(files_run, dataset + "/medium", "hdf5_results/Dvstar/distances_medium.hdf5", implementation, 8, 0)
        print("Large " + implementation + ".")
        res_our_large  = calculate_distance(files_run, dataset + "/large", "hdf5_results/Dvstar/distances_large.hdf5", implementation, 8, 0)
        
        compare_hdf5_files(dataset, "small", str(files_run))
        compare_hdf5_files(dataset, "medium", str(files_run))
        compare_hdf5_files(dataset, "large", str(files_run))
        
        catch_and_save(res_our_small, cwd / csv_filename, "small", files_run, th_small, min_small, max_small, implementation, 8)
        catch_and_save(res_our_medium, cwd / csv_filename, "medium", files_run, th_medium, min_medium, max_medium, implementation, 8)
        catch_and_save(res_our_large, cwd / csv_filename, "large", files_run, th_large, min_large, max_large, implementation, 8)
        
        files_run = files_run * 2 
        os.remove(cwd / "hdf5_results/Dvstar/distances_small.hdf5")
        os.remove(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")
        os.remove(cwd / "hdf5_results/Dvstar/distances_large.hdf5")

#####################################
# Function for benchmarking degree  #
# of parallelization on dataset.    #  
#####################################
def parallelization_benchmark(dataset: str, implementation: str):
    csv_filename = get_csv_name(dataset, "parallelization_")
    print(csv_filename)

    th_small, min_small, max_small = get_parameter_from_bintree(os.listdir(cwd / dataset / "small")[0])
    th_medium, min_medium, max_medium = get_parameter_from_bintree(os.listdir(cwd / dataset / "medium")[0])
    th_large, min_large, max_large = get_parameter_from_bintree(os.listdir(cwd / dataset / "large")[0])

    nb_files = count_nb_files(cwd / dataset / "small")

    files_run = 2
    while (files_run * 2 < nb_files):
        files_run = files_run * 2

    nr_cores = 1

    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_small.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_small.hdf5")
    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")
    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_large.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_large.hdf5")
    
    while (nr_cores <= 16):
        print("Benchmarking small with " + str(nr_cores) + " cores.")
        res_our_small  = calculate_distance(files_run, dataset + "/small", "hdf5_results/Dvstar/distances_small.hdf5", implementation, nr_cores, 0)
        print("Benchmarking medium with " + str(nr_cores) + " cores.")
        res_our_medium = calculate_distance(files_run, dataset + "/medium", "hdf5_results/Dvstar/distances_medium.hdf5", implementation, nr_cores, 0)
        print("Benchmarking large with " + str(nr_cores) + " cores.")
        res_our_large  = calculate_distance(files_run, dataset + "/large", "hdf5_results/Dvstar/distances_large.hdf5", implementation, nr_cores, 0)

        compare_hdf5_files(dataset, "small", str(files_run))
        compare_hdf5_files(dataset, "medium", str(files_run))
        compare_hdf5_files(dataset, "large", str(files_run))

        catch_and_save(res_our_small, cwd / csv_filename, "small", nb_files, th_small, min_small, max_small, implementation, nr_cores)
        catch_and_save(res_our_medium, cwd / csv_filename, "medium", nb_files, th_medium, min_medium, max_medium, implementation, nr_cores)
        catch_and_save(res_our_large, cwd / csv_filename, "large", nb_files, th_large, min_large, max_large, implementation, nr_cores)
        
        nr_cores = nr_cores * 2
        os.remove(cwd / "hdf5_results/Dvstar/distances_small.hdf5")
        os.remove(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")
        os.remove(cwd / "hdf5_results/Dvstar/distances_large.hdf5")

def compare_hdf5_files(dataset: str, size: str, files_run: str):
    vlmc_file = h5py.File(cwd / "hdf5_results" / "Dvstar" / ("distances_" + size + ".hdf5"), 'r')
    vlmc_group = vlmc_file.get("distances")
    vlmc_group = vlmc_group.get("distances")

    pst_file = h5py.File(cwd / "hdf5_results" / "PstClassifierSeqan" / (dataset + "/" + size + "_" + files_run + ".hdf5"), 'r')
    pst_group = pst_file.get("distances")
    pst_group = pst_group.get("distances")

    error_tolerance = 1e-7
    misses = 0

    for x in range(0, vlmc_group[0].size):
        for y in range(0, vlmc_group[0].size):
            if abs(vlmc_group[x][y] - pst_group[x][y]) > error_tolerance:
                misses = misses + 1
                print("Our = " + str(vlmc_group[x][y]) + " pst = " + str(pst_group[x][y]))
                print("Difference = " + str(vlmc_group[x][y] - pst_group[x][y]))

    if misses == 0:
        print("Correctness Check " + dataset + " : PASSED")

@app.command()
def benchmark():
    # ECOLI BENCHMARKING 
    # Pst_normal_benchmaking("data/benchmarking/ecoli")
    # normal_benchmaking("data/benchmarking/ecoli", "sorted-vector")
    # normal_benchmaking("data/benchmarking/ecoli", "hashmap")
    # normal_benchmaking("data/benchmarking/ecoli", "combo")

    # HUMAN BENCHMARKING
    # Pst_normal_benchmaking("data/benchmarking/human")
    # normal_benchmaking("data/benchmarking/human", "sorted-vector")
    # normal_benchmaking("data/benchmarking/human", "hashmap")
    # normal_benchmaking("data/benchmarking/human", "combo")

    parallelization_benchmark("data/benchmarking/human", "sorted-vector")

if __name__ == "__main__":
    app()
