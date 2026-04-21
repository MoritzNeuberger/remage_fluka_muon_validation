import argparse
from datetime import datetime
from pathlib import Path
import random
import shutil

import geometry as geom

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[1]
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "gen"


def build_input_file(
    filepath: Path,
    project_name: str,
    material: str,
    energy_gev: float,
    particle: str,
    n_primaries: int,
    burnin_z_cm: float = 400,
    seed: int = None,
) -> None:
    geom_data = geom.generate_geometry(material)

    beampos_z_cm = -geom_data["z_extent_cm"] / 2
    depth_cut_cm = burnin_z_cm
    z_entry_cm = -geom_data["z_extent_cm"] / 2
    z_exit_cm = geom_data["z_extent_cm"] / 2
    z_cut_cm = z_entry_cm + depth_cut_cm

    if seed is None:
        seed = random.randint(1, 999999)

    lines = [
        "FREE",
        f"* Auto-generated project: {project_name}",
        f"* Material={material}, particle={particle}, E={energy_gev:.6g} GeV, primaries={n_primaries}",
        "TITLE",
        f"Muon neutron multiplicity in {material} ({project_name})",
        "DEFAULTS 0.0 0.0 0.0 0.0 0.0 0.0 PRECISIO",
        f"BEAM {energy_gev:.6g} 0.0 0.0 0.0 0.0 0.0 {particle}",
        f"BEAMPOS 0.0 0.0 {beampos_z_cm:.6f}",
        "* Enable muon photonuclear and real-photon photonuclear on detector medium",
        "MUPHOTON 1.0 0.0 0.0 M000 M000 1.0",
        "PHOTONUC 1.0 0.0 0.0 M000 M000 1.0",
        "IONTRANS -2.0",
        "PHYSICS 3.0 0.0 0.0 0.0 0.0 0.0 EVAPORAT",
        geom_data["text"],
        "* User routine hooks for strict per-event neutron multiplicity",
        f"USRICALL {depth_cut_cm:.6f} {z_entry_cm:.6f} {z_exit_cm:.6f} {z_cut_cm:.6f} 0.0 0.0 NCOUNT",
        "USERDUMP 100.0 99.0 7.0 1.0 0.0 0.0 MGDRAW",
        "RESNUCLEI 3.0 0.0 0.0 0.0 -1",
        f"RANDOMIZE 1.0 {seed}",
        f"START {n_primaries}",
        "STOP",
    ]

    filepath.write_text("\n".join(lines), encoding="utf-8")


def parse_args(external_args: list = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a complete FLUKA project package for neutron and isotope production "
            "studies (input + companion user-routine template)."
        )
    )
    parser.add_argument(
        "--material",
        type=str,
        required=True,
        help="Material name for the FLUKA simulation.",
    )
    parser.add_argument(
        "--energy-gev",
        type=float,
        default=100.0,
        help="Primary muon kinetic energy in GeV (default: 100).",
    )
    parser.add_argument(
        "--beam-particle",
        type=str,
        default="MUON-",
        help="FLUKA BEAM particle name (default: MUON-).",
    )
    parser.add_argument(
        "--n-primaries",
        type=int,
        default=20000,
        help="Number of primary events on START card (default: 20000).",
    )
    parser.add_argument(
        "--depth-cut-cm",
        type=float,
        default=400.0,
        help="Depth cut from target entrance in centimeters (default: 400.0).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed value for RANDOMIZE card (default: 1).",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Root folder where project folder is created (default: gen/fluka_input).",
    )
    parser.add_argument(
        "--if-exists",
        choices=("fail", "overwrite", "timestamp"),
        default="overwrite",
        help="Behavior when project folder exists (default: fail).",
    )
    parser.add_argument(
        "--fluka-bin-path",
        type=Path,
        default=Path(REPO_ROOT / "sw" / "fluka4-5.1" / "bin"),
        help="Path to the FLUKA binary directory (default: sw/fluka4-5.1/bin).",
    )
    parser.add_argument(
        "--idx",
        type=int,
        default=None,
        help=(
            "Optional index to append to project folder name for uniqueness. "
            "If not provided, --if-exists behavior is applied instead."
        ),
    )
    return parser.parse_args(external_args)

def prepare_project_dir(output_root: Path, project_name: str, if_exists: str, idx: int) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    if idx is not None:
        project_dir = output_root / f"{project_name}_{idx}"
    else:
        project_dir = output_root / project_name
        
    if not project_dir.exists():
        project_dir.mkdir(parents=True)
        return project_dir

    if if_exists == "fail":
        msg = (
            f"Project folder already exists: {project_dir}. "
            "Use --if-exists overwrite or --if-exists timestamp."
        )
        raise FileExistsError(msg)

    if if_exists == "overwrite":
        shutil.rmtree(project_dir)
        project_dir.mkdir(parents=True)
        return project_dir

    suffix = datetime.now().strftime("%Y%m%d_%H%M%S")
    stamped_dir = output_root / f"{project_name}_{suffix}"
    stamped_dir.mkdir(parents=True)
    return stamped_dir

def generate_project(external_args: list = None):
    args = parse_args(external_args)

    if args.energy_gev <= 0:
        raise ValueError("--energy-gev must be > 0.")
    if args.n_primaries <= 0:
        raise ValueError("--n-primaries must be > 0.")
    if args.depth_cut_cm < 0:
        raise ValueError("--depth-cut-cm must be >= 0.")

    output_root = args.output_root.resolve()
    project_dir = prepare_project_dir(output_root, args.material, args.if_exists, args.idx)

    input_file_path = project_dir / "fluka_input.inp"
    build_input_file(
        filepath=input_file_path,
        project_name=args.material,
        material=args.material,
        energy_gev=args.energy_gev,
        particle=args.beam_particle,
        n_primaries=args.n_primaries,
        burnin_z_cm=args.depth_cut_cm,
        seed=args.seed,
    )
	
    fortran_file_path = project_dir / "neutron_scoring.f"
    shutil.copy(REPO_ROOT / "scripts" / "misc" / "neutron_scoring.f", fortran_file_path)
    
    compile_script_path = project_dir / "compile_and_run.sh"
    shutil.copy(REPO_ROOT / "scripts" / "misc" / "compile_and_run.sh", compile_script_path)
    with Path(compile_script_path).open("r", encoding="utf-8") as fh:
        compile_script_text = fh.read()
    compile_script_text = compile_script_text.replace("{rel_fluka_bin}", str(args.fluka_bin_path))
    with Path(compile_script_path).open("w", encoding="utf-8") as fh:
        fh.write(compile_script_text)
    return project_dir

if __name__ == "__main__":
    generate_project()