import project as proj
from pathlib import Path
import subprocess
from concurrent.futures import ProcessPoolExecutor

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[1]
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "gen"
N_THREADS = 20

def run_project(project_dir):    
	compile_script_path = project_dir / "compile_and_run.sh"
	print(f"Running FLUKA simulation in {project_dir}...")
	subprocess.run(["bash", str(compile_script_path)], cwd=project_dir, check=False)
    
def generate_projects(materials, n_projects_per_material=1, n_primaries=100000):
	project_dirs = {}
	for material in materials:
		project_dirs[material] = []
		for idx in range(n_projects_per_material):
			print(f"Generating project for {material}...")
			args = [
				"--material", material,
				"--energy-gev", "100",
				"--n-primaries", str(n_primaries),
				"--depth-cut-cm", "400",
				"--output-root", str(DEFAULT_OUTPUT_ROOT),
				"--if-exists", "overwrite",
			]
			if n_projects_per_material > 1:
				args += ["--idx", f"{idx}"]
			project_dir = proj.generate_project(args)
			project_dirs[material].append(project_dir)
	return project_dirs

def run_sims(materials, n_projects_per_material=1, n_threads=N_THREADS, n_primaries=100000):
	project_dirs = generate_projects(materials, n_projects_per_material=n_projects_per_material, n_primaries=n_primaries)
	project_paths = [project_dir for material_dirs in project_dirs.values() for project_dir in material_dirs]
	print(f"Running simulations for {len(project_paths)} projects with {n_threads} threads...")
	with ProcessPoolExecutor(max_workers=n_threads) as executor:
		executor.map(run_project, project_paths)
	return project_dirs	

