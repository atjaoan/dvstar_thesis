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
import random 

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
    
def get_csv_name(dataset_primary: str, dataset_secondary: str, now: datetime):
    out_path = "csv_results/" + now.strftime("%m_%d") + "/"
    #os.makedirs(cwd / out_path, exist_ok=True)
    return out_path + dataset_primary + "_" + dataset_secondary + "_" + now.strftime("%m_%d_%H_%M") + ".csv"

def count_nb_files(dir: str):
    count = 0
    for file in os.listdir(dir):
        count += 1
    return count 

##############################
# Get build parameters from  #
# size.                  #
##############################
def get_parameter_from_bintree(size: str) -> tuple[float, int, int]:
    if size == "small":
        return 3.9075, 9, 6
    elif size == "medium":
        return 3, 6, 8
    elif size == "large":
        return 2, 3, 10 
    else:
        return 3.9075, 10, 9

####################################
# Call to our implementation of    #
# distance calculation between     #
# one or two directories of VLMCs. #
####################################     
def calculate_distance(set_size: int, genome_path_primary: str, genome_path_secondary: str, out_path: str, vlmc_container: str, nr_cores: int, background_order: int) -> subprocess.CompletedProcess:
    args = (
        "perf", "stat",
        "-r",
        "3",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "build/dist", 
        "-p", cwd / genome_path_primary,
        "-s", cwd / genome_path_secondary,
        "-v", vlmc_container,
        "-n", str(nr_cores),
        "-b", str(background_order),
        "-a", str(set_size)# ,
        # "-o", str(out_path)
    )
    return subprocess.run(args, capture_output=True, text=True)

######################################
# Call to the old implementation of  #
# distance calculation between VLMCs #
# from PstClassifierSeqan submodule. #
######################################
def PstClassifierSeqan(set_size: int, genome_path_primary: str, genome_path_secondary: str, out_path: str, background_order: int) -> subprocess.CompletedProcess:
    args = (
        "perf", "stat",
        "-r",
        "3",
        "-e branch-misses,branches,task-clock,cycles,instructions,cache-references,cache-misses",
        cwd / "submodules/PstClassifierSeqan/build/src/calculate-distances", 
        "-p", cwd / genome_path_primary,
        "-y", cwd / genome_path_secondary,
        "-n", "dvstar",
        "-a", str(set_size),
        "-b", str(background_order)# ,
        # "-s", cwd / out_path
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

    data = [get_git_commit_version(), implementation, vlmc_size, set_size, threshold, min_count, max_depth, nr_cores_used, combo_init_size]
    columns = ["repo_version", "implementation", "vlmc_size", "set_size", "threshold", "min_count", "max_depth", "nr_cores_used", "combo_init_size"]
    
    for line in new_line_separated_attr:
        split_line = line.split('#')

        divergence = split_line[1].split('(')[1].split('-')[1].split('%')[0]
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

        # try: 
        #     data.extend([float(right_value), int(float(count)), float(divergence)])
        # except:
        #     print("right_value = " + right_value)
        #     print("count = " + count)
        #     print("divergence = " + divergence)
        #     raise Exception("Failed to write to csv")
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
def Pst_normal_benchmaking(path_primary: str, path_secondary: str, csv_filename: str, single_run: bool, set_size: int):
    dset_p = path_primary.split("/")
    dataset_primary = dset_p[len(dset_p) - 1]

    dset_s = path_secondary.split("/")
    dataset_secondary = dset_s[len(dset_s) - 1]

    nb_files = min(count_nb_files(cwd / path_primary / "small"), count_nb_files(cwd / path_secondary / "small"))

    files_run = 32  

    th_small, min_small, max_small = get_parameter_from_bintree("small")
    th_medium, min_medium, max_medium = get_parameter_from_bintree("medium")
    th_large, min_large, max_large = get_parameter_from_bintree("large")

    out_path = "hdf5_results/PstClassifierSeqan/" + dataset_primary + "_to_" + dataset_secondary + "/"
    
    #os.makedirs(out_path, exist_ok=True)

    if single_run:
        print("Small PstClassifierSeqan.")
        res_small  = PstClassifierSeqan(set_size, path_primary + "/small", path_secondary + "/small", out_path + "small_" + str(files_run) + ".hdf5", 0)
        print(res_small.stderr)
        print("Medium PstClassifierSeqan.")
        res_medium = PstClassifierSeqan(set_size, path_primary + "/medium", path_secondary + "/medium", out_path + "medium_" + str(files_run) + ".hdf5", 0)
        print(res_medium.stderr)
        print("Large PstClassifierSeqan.")
        res_large  = PstClassifierSeqan(set_size, path_primary + "/large", path_secondary + "/large", out_path + "large_" + str(files_run) + ".hdf5", 0)
        print(res_large.stderr)

        catch_and_save(res_small, cwd / csv_filename, "small", files_run, th_small, min_small, max_small, "PstClassifierSeqan", 8)
        catch_and_save(res_medium, cwd / csv_filename, "medium", files_run, th_medium, min_medium, max_medium, "PstClassifierSeqan", 8)
        catch_and_save(res_large, cwd / csv_filename, "large", files_run, th_large, min_large, max_large, "PstClassifierSeqan", 8)
    else:
        while((files_run < 10000) & (files_run <= nb_files)):
            print("Benchmarking with " + str(files_run) + " VLMCs...")
            print("Small PstClassifierSeqan.")
            res_small  = PstClassifierSeqan(files_run, path_primary + "/small", path_secondary + "/small", out_path + "small_" + str(files_run) + ".hdf5", 0)
            print("Medium PstClassifierSeqan.")
            res_medium = PstClassifierSeqan(files_run, path_primary + "/medium", path_secondary + "/medium", out_path + "medium_" + str(files_run) + ".hdf5", 0)
            print("Large PstClassifierSeqan.")
            res_large  = PstClassifierSeqan(files_run, path_primary + "/large", path_secondary + "/large", out_path + "large_" + str(files_run) + ".hdf5", 0)

            catch_and_save(res_small, cwd / csv_filename, "small", files_run, th_small, min_small, max_small, "PstClassifierSeqan", 8)
            catch_and_save(res_medium, cwd / csv_filename, "medium", files_run, th_medium, min_medium, max_medium, "PstClassifierSeqan", 8)
            catch_and_save(res_large, cwd / csv_filename, "large", files_run, th_large, min_large, max_large, "PstClassifierSeqan", 8)
            files_run = files_run * 2

######################################
# Function for benchmarking time and # 
# cache-misses on dataset with new   #
# implementation.                    #
######################################
def normal_benchmarking(path_primary: str, path_secondary: str, implementation: str, csv_filename: str, single_run: bool, cores: int, set_size: int):
    dset_p = path_primary.split("/")
    dataset_primary = dset_p[len(dset_p) - 1]

    dset_s = path_secondary.split("/")
    dataset_secondary = dset_s[len(dset_s) - 1]

    nb_files = min(count_nb_files(cwd / path_primary / "small"), count_nb_files(cwd / path_secondary / "small"))

    th_small, min_small, max_small = get_parameter_from_bintree("small")
    th_medium, min_medium, max_medium = get_parameter_from_bintree("medium")
    th_large, min_large, max_large = get_parameter_from_bintree("large")

    files_run = 32

    #out_path = "hdf5_results/Dvstar/"

    #os.makedirs(out_path, exist_ok=True)

    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_small.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_small.hdf5")
    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")
    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_large.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_large.hdf5")

    if single_run:
        print("Small " + implementation + ".")
        res_our_small  = calculate_distance(set_size, path_primary + "/small", path_secondary + "/small", "hdf5_results/Dvstar/distances_small.hdf5", implementation, cores, 0)
        print(res_our_small.stderr)
        print("Medium " + implementation + ".")
        res_our_medium = calculate_distance(set_size, path_primary + "/medium", path_secondary + "/medium", "hdf5_results/Dvstar/distances_medium.hdf5", implementation, cores, 0)
        print(res_our_medium.stderr)
        print("Large " + implementation + ".")
        res_our_large  = calculate_distance(set_size, path_primary + "/large", path_secondary + "/large", "hdf5_results/Dvstar/distances_large.hdf5", implementation, cores, 0)
        print(res_our_large.stderr)
        
        catch_and_save(res_our_small, cwd / csv_filename, "small", files_run, th_small, min_small, max_small, implementation, cores)
        catch_and_save(res_our_medium, cwd / csv_filename, "medium", files_run, th_medium, min_medium, max_medium, implementation, cores)
        catch_and_save(res_our_large, cwd / csv_filename, "large", files_run, th_large, min_large, max_large, implementation, cores)
    else:
        while((files_run < 10000) & (files_run <= nb_files)):
            print("Benchmarking with " + str(files_run) + " VLMCs...")

            print("Small " + implementation + ".")
            res_our_small  = calculate_distance(files_run, path_primary + "/small", path_secondary + "/small", "hdf5_results/Dvstar/distances_small.hdf5", implementation, cores, 0)
            print("Medium " + implementation + ".")
            res_our_medium = calculate_distance(files_run, path_primary + "/medium", path_secondary + "/medium", "hdf5_results/Dvstar/distances_medium.hdf5", implementation, cores, 0)
            print("Large " + implementation + ".")
            res_our_large  = calculate_distance(files_run, path_primary + "/large", path_secondary + "/large", "hdf5_results/Dvstar/distances_large.hdf5", implementation, cores, 0)

            #compare_hdf5_files(dataset_primary + "_to_" + dataset_secondary, "small", str(files_run))
            #compare_hdf5_files(dataset_primary + "_to_" + dataset_secondary, "medium", str(files_run))
            #compare_hdf5_files(dataset_primary + "_to_" + dataset_secondary, "large", str(files_run))

            catch_and_save(res_our_small, cwd / csv_filename, "small", files_run, th_small, min_small, max_small, implementation, cores)
            catch_and_save(res_our_medium, cwd / csv_filename, "medium", files_run, th_medium, min_medium, max_medium, implementation, cores)
            catch_and_save(res_our_large, cwd / csv_filename, "large", files_run, th_large, min_large, max_large, implementation, cores)

            files_run = files_run * 2 
            #os.remove(cwd / "hdf5_results/Dvstar/distances_small.hdf5")
            #os.remove(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")
            #os.remove(cwd / "hdf5_results/Dvstar/distances_large.hdf5")

#####################################
# Function for benchmarking degree  #
# of parallelization on dataset.    #  
#####################################
def parallelization_benchmark(dataset: str, implementation: str, max_cores: int, csv_filename: str, set_size: int = -1):
    th_small, min_small, max_small = get_parameter_from_bintree("small")
    th_medium, min_medium, max_medium = get_parameter_from_bintree("medium")
    th_large, min_large, max_large = get_parameter_from_bintree("large")

    nb_files = count_nb_files(cwd / dataset / "small")

    nr_cores = 1

    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_small.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_small.hdf5")
    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")
    if (os.path.isfile(cwd / "hdf5_results/Dvstar/distances_large.hdf5")):
        os.remove(cwd / "hdf5_results/Dvstar/distances_large.hdf5")
    
    while (nr_cores <= max_cores):
        print("Benchmarking small with " + str(nr_cores) + " cores.")
        res_our_small  = calculate_distance(set_size, dataset + "/small", dataset + "/small", "hdf5_results/Dvstar/distances_small.hdf5", implementation, nr_cores, 0)
        print("Benchmarking medium with " + str(nr_cores) + " cores.")
        res_our_medium = calculate_distance(set_size, dataset + "/medium", dataset + "/medium", "hdf5_results/Dvstar/distances_medium.hdf5", implementation, nr_cores, 0)
        print("Benchmarking large with " + str(nr_cores) + " cores.")
        res_our_large  = calculate_distance(set_size, dataset + "/large", dataset + "/large", "hdf5_results/Dvstar/distances_large.hdf5", implementation, nr_cores, 0)

        #compare_hdf5_files(dataset, "small", str(files_run))
        #compare_hdf5_files(dataset, "medium", str(files_run))
        #compare_hdf5_files(dataset, "large", str(files_run))

        catch_and_save(res_our_small, cwd / csv_filename, "small", nb_files, th_small, min_small, max_small, implementation, nr_cores)
        catch_and_save(res_our_medium, cwd / csv_filename, "medium", nb_files, th_medium, min_medium, max_medium, implementation, nr_cores)
        catch_and_save(res_our_large, cwd / csv_filename, "large", nb_files, th_large, min_large, max_large, implementation, nr_cores)
        
        if nr_cores == 1:
            nr_cores = 2
        else:
            nr_cores += 2
        #os.remove(cwd / "hdf5_results/Dvstar/distances_small.hdf5")
        #os.remove(cwd / "hdf5_results/Dvstar/distances_medium.hdf5")
        #os.remove(cwd / "hdf5_results/Dvstar/distances_large.hdf5")

def compare_hdf5_files(dataset: str, size: str, files_run: str):
    vlmc_file = h5py.File(cwd / "hdf5_results" / "Dvstar" / ("distances_" + size + ".hdf5"), 'r')
    vlmc_group = vlmc_file.get("distances")
    vlmc_group = vlmc_group.get("distances")

    pst_file = h5py.File(cwd / "hdf5_results" / "PstClassifierSeqan" / (dataset + "/" + size + "_" + files_run + ".hdf5"), 'r')
    pst_group = pst_file.get("distances")
    pst_group = pst_group.get("distances")

    error_tolerance = 1e-7
    misses = 0

    for x, row in enumerate(pst_group):
        for y, _ in enumerate(row):
            if abs(vlmc_group[x][y] - pst_group[x][y]) > error_tolerance:
                misses = misses + 1
                print("Our = " + str(vlmc_group[x][y]) + " pst = " + str(pst_group[x][y]))
                print("Difference = " + str(vlmc_group[x][y] - pst_group[x][y]))

    if misses != 0:
        raise Exception("Correctness Check " + dataset + " with " + size + " size and " + files_run + " files : FAILED")

@app.command()
def benchmark(single_run: bool = True, cores: int = 8, set_size_virus: int = -1):
    now = datetime.now()
    datasets = [("human", "human"), ("human", "turkey"), ("human", "corn"), ("turkey", "turkey"), ("turkey", "corn"), ("corn", "corn"), ("ecoli", "ecoli")]
    containers = ["sorted-vector", "sorted-search", "hashmap", "veb", "ey", "alt-btree"] # ,"sorted-search-ey"]
    for primary, secondary in datasets:
        print("Benchmarking " + primary + " to " + secondary)
        csv_filename = get_csv_name(primary, secondary, now)
        p = "data/benchmarking/" + primary
        s = "data/benchmarking/" + secondary
        if primary == "ecoli" or secondary == "ecoli":
            Pst_normal_benchmaking(p, s, csv_filename, single_run, set_size_virus)
        else:
            Pst_normal_benchmaking(p, s, csv_filename, single_run, -1)
        for container in containers:
            if primary == "ecoli" or secondary == "ecoli":
                normal_benchmarking(p, s, container, csv_filename, single_run, cores, set_size_virus)
            else :
                normal_benchmarking(p, s, container, csv_filename, single_run, cores, -1)

@app.command()
def sizebenchmark(single_run: bool = False, cores: int = 8):
    now = datetime.now()
    datasets = [("ecoli", "ecoli")]
    containers = ["sorted-vector", "sorted-search", "hashmap", "veb", "ey", "alt-btree"] # ,"sorted-search-ey"]
    for primary, secondary in datasets:
        print("Benchmarking " + primary + " to " + secondary)
        csv_filename = get_csv_name(primary, secondary, now)
        p = "data/benchmarking/" + primary
        s = "data/benchmarking/" + secondary
        Pst_normal_benchmaking(p, s, csv_filename, single_run, -1)
        for container in containers:
            normal_benchmarking(p, s, container, csv_filename, single_run, cores, -1)

@app.command()
def parabench(cores: int = 8, species: str = "human", set_size: int = -1):
    containers = ["sorted-vector", "sorted-search", "hashmap", "veb", "ey", "alt-btree"]
    now = datetime.now()
    dataset = "data/benchmarking/" + species + "/"
    csv_filename = get_csv_name(dataset.split("/")[-2], "parallelization", now)
    for container in containers:
        parallelization_benchmark(dataset, container, cores, csv_filename, set_size)

# @app.command()
# def change_names():
#     dir = cwd / "data/benchmarking/bigdata/large" 
#     list_paths = []
#     for path in os.listdir(cwd / "data/benchmarking/bigdata/large"):
#         list_paths.append(path)
# 
#     new_names = random.sample(range(1, 1000), len(list_paths))
#     for i, path in enumerate(list_paths): 
#         os.rename(dir / path, dir / (str(new_names[i]) + ".bintree"))

if __name__ == "__main__":
    app()
