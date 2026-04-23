import project as proj
from pathlib import Path
import subprocess
from concurrent.futures import ProcessPoolExecutor
import json

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[1]
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "gen"
N_THREADS = 20


def _looks_like_project_dir(project_dir: Path) -> bool:
	if not project_dir.is_dir():
		return False

	# Dual-mode layout: project root contains neutron/ and isotope/ mode subfolders.
	if (project_dir / "neutron").is_dir() or (project_dir / "isotope").is_dir():
		return True

	# Legacy single-folder layout fallback.
	if (project_dir / "fluka_input.inp").exists():
		return True

	return False


def discover_existing_projects(cases, output_root=DEFAULT_OUTPUT_ROOT):
	output_root = Path(output_root)
	project_dirs = {}

	if not output_root.exists():
		print(f"[warn] output root does not exist: {output_root}")
		return project_dirs

	for case in cases:
		material, energy = case
		base_name = f"{material}-{energy}GeV"
		case_dirs = []

		base_dir = output_root / base_name
		if _looks_like_project_dir(base_dir):
			case_dirs.append(base_dir)

		indexed_dirs = []
		prefix = f"{base_name}_"
		for candidate in output_root.iterdir():
			if not candidate.name.startswith(prefix):
				continue
			suffix = candidate.name[len(prefix):]
			if not suffix.isdigit():
				continue
			if not _looks_like_project_dir(candidate):
				continue
			indexed_dirs.append((int(suffix), candidate))

		indexed_dirs.sort(key=lambda item: item[0])
		case_dirs.extend(path for _, path in indexed_dirs)

		if not case_dirs:
			print(f"[warn] no existing project folders found for case: {material}, {energy} GeV")
			continue

		project_dirs[case] = case_dirs
		print(f"Discovered {len(case_dirs)} existing project folders for case: {material}, {energy} GeV")

	return project_dirs

def _run_mode(project_dir: Path, mode: str):
	mode_dir = project_dir / mode
	if not mode_dir.exists():
		return {"status": "missing", "returncode": None, "path": str(mode_dir)}

	compile_script_path = mode_dir / "compile_and_run.sh"
	if not compile_script_path.exists():
		return {"status": "missing-script", "returncode": None, "path": str(mode_dir)}

	print(f"Running FLUKA simulation in {mode_dir}...")
	completed = subprocess.run(["bash", str(compile_script_path)], cwd=mode_dir, check=False)
	status = "success" if completed.returncode == 0 else "failed"
	return {"status": status, "returncode": completed.returncode, "path": str(mode_dir)}


def run_project(project_dir):
	status = {
		"project": str(project_dir),
		"modes": {
			"neutron": _run_mode(project_dir, "neutron"),
			"isotope": _run_mode(project_dir, "isotope"),
		},
	}

	status_file = project_dir / "run_status.json"
	with status_file.open("w", encoding="utf-8") as fh:
		json.dump(status, fh, indent=2, sort_keys=True)
		fh.write("\n")
	return status
    
def generate_projects(cases, n_projects_per_case=1, n_primaries=100000):
	project_dirs = {}
	for case in cases:
		material = case[0]
		energy = case[1]
		project_dirs[case] = []
		for idx in range(n_projects_per_case):
			print(f"Generating dual-mode project for {material}...")
			seed = int((energy * 10000) + (idx + 1) + (sum(ord(c) for c in material) * 13))
			args = [
				"--material", material,
				"--energy-gev", str(energy),
				"--n-primaries", str(n_primaries),
				"--depth-cut-cm", "400",
				"--physics-mode", "both",
				"--seed", str(seed),
				"--output-root", str(DEFAULT_OUTPUT_ROOT),
				"--if-exists", "overwrite",
				"--name", f"{material}-{energy}GeV"
			]
			if n_projects_per_case > 1:
				args += ["--idx", f"{idx}"]
			project_dir = proj.generate_project(args)
			project_dirs[case].append(project_dir)
	return project_dirs

def run_sims(cases, n_projects_per_case=1, n_threads=N_THREADS, n_primaries=100000):
	project_dirs = generate_projects(cases, n_projects_per_case=n_projects_per_case, n_primaries=n_primaries)
	project_paths = [project_dir for case_dirs in project_dirs.values() for project_dir in case_dirs]
	print(f"Running simulations for {len(project_paths)} projects with {n_threads} threads...")
	with ProcessPoolExecutor(max_workers=n_threads) as executor:
		results = list(executor.map(run_project, project_paths))

	failed_modes = 0
	for result in results:
		for mode_result in result["modes"].values():
			if mode_result["status"] == "failed":
				failed_modes += 1

	if failed_modes > 0:
		print(f"Completed with {failed_modes} failed mode runs. See run_status.json files in each project folder.")
	else:
		print("Completed all mode runs successfully.")
	return project_dirs	

