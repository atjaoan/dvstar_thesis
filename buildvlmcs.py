import typer
import subprocess
from pathlib import Path
import os 

app = typer.Typer()

cwd = Path(__file__).parent

def get_bintree_name(genome_path: str, threshold: float, min_count: int, max_depth: int):
    return os.path.splitext(genome_path)[0] + f"_{threshold}_{min_count}_{max_depth}.bintree"

############################################
# Builds VLMCs from sequences.             #
# Dvstars projects implementation is used. #
############################################ 
def dvstar_build(genome_path: Path, out_path: Path, threshold: float, min_count: int, max_depth: int):
    for genome in os.listdir(genome_path):
        args = (
            "./build/dvstar",
            "--mode", "build",
            "--threshold", str(threshold),
            "--min-count", str(min_count),
            "--max-depth", str(max_depth),
            "--fasta-path", genome_path / genome,
            "--out-path", out_path / get_bintree_name(genome, threshold, min_count, max_depth)
        )
        subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

@app.command()
def build():
    print("Building")

if __name__ == "__main__":
    app()