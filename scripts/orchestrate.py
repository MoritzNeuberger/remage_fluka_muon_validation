import run_sims as sims
import post_proc as postp
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Orchestrate FLUKA simulations and post-processing.")
    parser.add_argument("--n-projects-per-material", type=int, default=20, help="Number of projects to generate per material (default: 20).")
    parser.add_argument("--n-threads", type=int, default=20, help="Number of parallel threads for running simulations (default: 20).")
    parser.add_argument("--n-primaries", type=int, default=100000, help="Number of primary events for each simulation (default: 100000).")
    return parser.parse_args()

materials = ["lar", "water", "rock"]

if __name__ == "__main__":
    args = parse_args()
    project_folders = sims.run_sims(materials, n_projects_per_material=args.n_projects_per_material, n_threads=args.n_threads, n_primaries=args.n_primaries )
    postp.post_proc(project_folders)