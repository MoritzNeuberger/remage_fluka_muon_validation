import run_sims as sims
import post_proc as postp
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Orchestrate FLUKA simulations and post-processing.")
    parser.add_argument("--n-projects-per-case", type=int, default=20, help="Number of projects to generate per case (default: 20).")
    parser.add_argument("--n-threads", type=int, default=20, help="Number of parallel threads for running simulations (default: 20).")
    parser.add_argument("--n-primaries", type=int, default=100000, help="Number of primary events for each simulation (default: 100000).")
    parser.add_argument(
        "--post-proc-only",
        action="store_true",
        help="Skip simulations and only run post-processing from existing project folders in gen/.",
    )
    parser.add_argument("--specific-case", type=str, help="Run only a specific case in the format 'material-energy' (e.g., 'lar-100'). Overrides --n-projects-per-case.")
    return parser.parse_args()

cases = [
    ("lar", 100),
    ("lar", 280),
    ("water", 100),
    ("water", 280),
    ("rock", 100),
    ("enrGe", 100)
]

if __name__ == "__main__":
    args = parse_args()
    if args.specific_case:
        material, energy = args.specific_case.split("-")
        cases = [(material, int(energy))]
    if args.post_proc_only:
        print("Post-processing only mode enabled.")
        print("Ignoring --n-projects-per-case, --n-threads, and --n-primaries.")
        project_folders = sims.discover_existing_projects(cases)
        n_discovered = sum(len(folders) for folders in project_folders.values())
        if n_discovered == 0:
            raise SystemExit("No existing project folders were found for configured cases. Nothing to post-process.")
    else:
        project_folders = sims.run_sims(
            cases,
            n_projects_per_case=args.n_projects_per_case,
            n_threads=args.n_threads,
            n_primaries=args.n_primaries,
        )
    postp.post_proc(project_folders)