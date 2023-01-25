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

def build_vlmc(kmc_db: Path) -> Path:
    args = (
        "build/build_vlmc",
        "--mode",
        "build-from-kmc-db",
        "--in-path",
        kmc_db,
        "--out-path",
        kmc_db.with_suffix(".bintree"),
        "--max-depth",
        f"{K}",
    )
    subprocess.run(args)

    return kmc_db.with_suffix(".bintree")


def count_kmers(fasta_path: Path) -> Path:
    args = (
        cwd / "build/kmc",
        "-b",
        "-ci1",
        "-cs4294967295",
        "-r",
        "-m24",
        f"-k{K + 1}",
        "-fm",
        fasta_path,
        f"tmp/{fasta_path.stem}",
        "tmp/",
    )
    subprocess.run(args)

    return Path("tmp") / fasta_path.stem

def dvstar(fasta_path: Path, out_path: Path) -> tuple[Path, float, float]:
    args = (
        cwd / "build/dvstar",
        "--fasta-path",
        fasta_path,
        "--threshold",
        "3.9075",
        "--max-depth",
        "4",
        "--min-count",
        "100",
        "--out-path",
        out_path, # <- change this
        "--temp-path",
        "tmp"
    )
    subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return out_path

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

def save_to_csv(res: subprocess.CompletedProcess, csv_path: Path):
    new_line_separated_attr = res.stderr.split('\n')[3:-8]

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

@app.command()
def run():
    fasta_path_1 = cwd / "tests/NC_001497.2.fa"
    fasta_path_2 = cwd / "tests/NC_028367.1.fa"
    vlmc_1 = dvstar(fasta_path_1, Path("NC_022098.1.bintree"))
    vlmc_2 = dvstar(fasta_path_2,  Path("NC_022099.1.bintree"))

    res = dvstar_cmp(vlmc_1, vlmc_2)

    save_to_csv(res, cwd / "tmp/benchmarks/test.csv")

if __name__ == "__main__":
    app()