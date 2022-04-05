import subprocess
import time
from pathlib import Path
from typing import Final

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import typer

app = typer.Typer()

K: Final[int] = 9


def convert_to_fasta(path: Path) -> Path:
    if path.suffix == ".fastq":
        if not path.with_suffix(".fasta").is_file():
            reference_convert_args = f"./seqtk seq -a {path}"
            r = subprocess.run(reference_convert_args.split(" "), capture_output=True)
            path = path.with_suffix(".fasta")
            with path.open("w") as f:
                f.write(r.stdout.decode())
        else:
            path = path.with_suffix(".fasta")
    return path


def grow_and_build(fasta_path: Path, percentage_to_add: float = 0.1) -> tuple[Path, float, float]:
    fasta_path = convert_to_fasta(fasta_path)

    records = SeqIO.parse(fasta_path, "fasta")
    seq = "".join(str(rec.seq) for rec in records)
    short_seq = seq[:int(len(seq) * percentage_to_add)]
    new_record = SeqRecord(id=fasta_path.stem + "-larger", seq=Seq(seq + short_seq))

    tmp_path = Path("tmp") / f"grown_{fasta_path.stem}_{int(percentage_to_add * 100)}_percent.fasta"
    with tmp_path.open("w") as f:
        f.write(new_record.format("fasta"))

    start = time.perf_counter()
    kmc_db, _, _ = count_kmers(tmp_path)
    path, _, _ = build_vlmc(kmc_db)
    end = time.perf_counter()
    return path, start, end


def add_10_percent(fasta_path: Path, kmc_db: Path) -> tuple[Path, float, float]:
    fasta_path = convert_to_fasta(fasta_path)

    records = SeqIO.parse(fasta_path, "fasta")
    seq = "".join(str(rec.seq) for rec in records)
    short_seq = seq[:len(seq) // 10]
    new_record = SeqRecord(id="Short fasta", seq=Seq(short_seq))

    tmp_path = Path("tmp") / "plus_10_percent.fasta"
    with tmp_path.open("w") as f:
        f.write(new_record.format("fasta"))

    start = time.perf_counter()
    kmc_db_new, _, _ = count_kmers(tmp_path)
    path, union_start, union_end = join_kmc_dbs(kmc_db, kmc_db_new)
    end = time.perf_counter()

    print(f"\tUnion takes: {union_end - union_start}")

    return path, start, end


def join_kmc_dbs(kmc_db1: Path, kmc_db2: Path) -> tuple[Path, float, float]:
    new_db_path = kmc_db1.parent / f"{kmc_db1.stem}_{kmc_db2.stem}_union"
    args = (
        "build/kmc_tools",
        "simple",
        kmc_db1,
        kmc_db2,
        "union",
        new_db_path,
        "-ocsum",
    )
    start = time.perf_counter()
    subprocess.run(args)
    end = time.perf_counter()
    print(" ".join([str(v) for v in args]))

    return new_db_path, start, end


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


@app.command()
def run(fasta_path: Path):
    kmc_db, start_kmc, stop_kmc = count_kmers(fasta_path)

    _, start_vlmc, stop_vlmc = build_vlmc(kmc_db)

    new_db, start_add_kmc, stop_add_kmc = add_10_percent(fasta_path, kmc_db)

    _, start_union_vlmc, stop_union_vlmc = build_vlmc(new_db)

    _, start_grow, stop_grow = grow_and_build(fasta_path, 0.1)


    print(f"For {fasta_path.stem}, max depth = {K}, file size = {fasta_path.stat().st_size}")
    print(f"\tKmc takes:\t\t {stop_kmc - start_kmc}s")
    print(f"\tVLMC takes:\t\t {stop_vlmc - start_vlmc}s")
    print(f"\tAdding 10 percent takes: {stop_add_kmc - start_add_kmc}s")
    print(f"\tBuilding on union takes: {stop_union_vlmc - start_union_vlmc}s")

    print(f"\n\tTotal original construction:\t\t {stop_vlmc - start_kmc}")
    print(f"\tTotal add 10 percent:\t\t\t {stop_union_vlmc - start_add_kmc}")
    print(f"\tBuilding on 10 percent larger file:\t {stop_grow - start_grow}s")


@app.command()
def run_default():
    pandoravirus_path = Path("python-prototype/NC_022098.1.fasta")
    downloads_folder = Path.home() / "data" / "ecoli"
    aid = "NC_002695.2.fasta"
    ecoli_path = downloads_folder / aid

    for fasta_path in [pandoravirus_path, ecoli_path]:
        run(fasta_path)


if __name__ == "__main__":
    app()
