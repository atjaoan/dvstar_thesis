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

def build_vlmc(kmc_db: Path) -> tuple[Path, float, float]:
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
    start = time.perf_counter()
    subprocess.run(args)
    end = time.perf_counter()

    return kmc_db.with_suffix(".bintree"), start, end


def count_kmers(fasta_path: Path) -> tuple[Path, float, float]:
    args = (
        "build/kmc",
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
    start = time.perf_counter()
    subprocess.run(args)
    end = time.perf_counter()

    return Path("tmp") / fasta_path.stem, start, end

def dvstar(fasta_path: Path, out_path: Path) -> tuple[Path, float, float]:
    args = (
        "build/dvstar",
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
        "build/dvstar",
        "--mode",
        "5",
        "--in-path",
        path_1,
        "--to-path",
        path_2
    )
    res = subprocess.run(args, capture_output=True, text=True)

    return res

def save_to_txt(res: subprocess.CompletedProcess, csv_path: Path):
    split = res.stderr.split('\n')

    data = []
    columns = []
    
    for s in split:
        x = s.split('#')
        if len(x) > 1:
            val = x[1].lstrip().split(' ')[0]
            x = x[0].strip()
            if "msec" in x:
                y = x.split('msec')
            else: 
                y = re.split(r"\s{2,}", x)

            col  = y[1].replace(":u", "")
            cnts = y[0]

            data.extend([val, cnts])
            columns.extend([col, col + "_counts"])

    if not os.path.exists(csv_path):
        df = pd.DataFrame(columns=columns)
    else:
        df = pd.read_csv(csv_path, dtype=str)

    df.loc[len(df)] = data
    df.to_csv(csv_path, index=False)

@app.command()
def run(): ## fasta_path: Path):
    fasta_path_1 = Path("/home/holmse/thesis/dvstar_thesis/tests/NC_001497.2.fa")
    fasta_path_2 = Path("/home/holmse/thesis/dvstar_thesis/tests/NC_028367.1.fa")
    vlmc_1 = dvstar(fasta_path_1, Path("NC_022098.1.bintree"))
    vlmc_2 = dvstar(fasta_path_2,  Path("NC_022099.1.bintree"))

    res = dvstar_cmp(vlmc_1, vlmc_2)

    save_to_txt(res, "tmp/benchmarks/test.csv")

if __name__ == "__main__":
    app()