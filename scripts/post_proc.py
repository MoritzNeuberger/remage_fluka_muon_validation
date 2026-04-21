import json
import awkward as ak
import numpy as np
from dataclasses import dataclass
import re
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[1]
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "gen"
NEUTRON_SCORING_FILE = "fluka_input001_neutron_scoring.dat"
FLUKA_INPUT_FILE = "fluka_input.inp"
RESNUCLEI_OUTPUT_FILE = "fluka_input001.out"

DENSITIES = {
    "lar": 1.396,  # g/cm^3
    "water": 1.0,  # g/cm^3
    "rock": 2.65,  # g/cm^3
}
ELEMENTS_PER_MATERIAL = {
    "lar": [("Ar", 18)],
    "water": [("H", 1), ("O", 8)],
    "rock": [("O", 8), ("Si", 14)],  # O, Si
}

############# Neutron Yield Calculation #############

def get_neutron_production_data_per_folder(folder):
    folder = Path(folder)
    with (folder / NEUTRON_SCORING_FILE).open() as f:
        lines = f.readlines()
    data = []
    for line in lines:        
        if line.startswith("#"):
            continue
        line_split = line.strip().split()
        try:
            tmp = {
                    "evtID": int(line_split[0]),
                    "n_accepted": int(line_split[1]),
                    "n_total": int(line_split[2])
                }
            data.append(tmp)
        except (ValueError, IndexError):
            continue
    return ak.Array(data)

def get_target_length(folder):
    folder = Path(folder)
    with (folder / FLUKA_INPUT_FILE).open() as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("USRICALL"):
            xmin = float(line.strip().split()[3])
            xcut = float(line.strip().split()[4])
            xmax = float(line.strip().split()[2])
            return {
                "length_accepted": xmin - xcut, # cm
                "length_total": xmin - xmax, # cm
                "xmin": xmin, # cm
                "xcut": xcut, # cm
                "xmax": xmax # cm
            }

def calculate_neutron_yield(n_neutrons, length, n_muons, material):
    density = DENSITIES[material]
    return n_neutrons / (length * density * n_muons)

def calculate_neutron_yield_per_folder(folder):
    target_length = get_target_length(folder)
    neutron_production_data = get_neutron_production_data_per_folder(folder)
    n_muons = len(neutron_production_data)
    material = Path(folder).name.split("_")[0]
    output = {}
    output["yield_accepted"] = calculate_neutron_yield(
        np.sum(neutron_production_data["n_accepted"]),
        target_length["length_accepted"],
        n_muons,
        material
    )
    output["yield_total"] = calculate_neutron_yield(
        np.sum(neutron_production_data["n_total"]),
        target_length["length_total"],
        n_muons,
        material
    )
    return output

def calculate_neutron_yield_per_project(folders):
    yields = [calculate_neutron_yield_per_folder(folder) for folder in folders]
    yield_accepted = {"val":np.mean([y["yield_accepted"] for y in yields]),"unc":np.std([y["yield_accepted"] for y in yields])/np.sqrt(len(yields))}
    yield_total = {"val":np.mean([y["yield_total"] for y in yields]),"unc":np.std([y["yield_total"] for y in yields])/np.sqrt(len(yields))}
    return {
        "yield_accepted": yield_accepted,
        "yield_total": yield_total,
        "raw": yields
    }

def process_neutron_yield_material(project_folders, output_file=None):
    project_yields = calculate_neutron_yield_per_project(project_folders)
    if output_file:
        with open(output_file, "w") as f:
            json.dump(project_yields, f)
    return project_yields


############# Isotope Production #############

@dataclass
class RESNCULEi_output:
    z_axis: np.ndarray
    n_axis: np.ndarray
    values: np.ndarray

def parse_resnuclei_output(file):
    lines = Path(file).read_text(errors="ignore").splitlines()

    start_index = None
    end_index = None
    max_z = None
    max_nz = None
    min_nz = None
    k = None

    # Find header and matrix start
    for i, line in enumerate(lines):
        if "Data follow in a matrix" in line:
            start_index = i + 1

            # Parse k from: A(z,n-z-k), k: -5
            km = re.search(r"k:\s*(-?\d+)", line)
            if km:
                k = int(km.group(1))

            # Previous line has Max Z / Max N-Z / Min N-Z
            pm = re.search(
                r"Max\.\s*Z:\s*(\d+),\s*Max\.\s*N-Z:\s*(-?\d+)\s*Min\.\s*N-Z:\s*(-?\d+)",
                lines[i - 1],
            )
            if pm:
                max_z = int(pm.group(1))
                max_nz = int(pm.group(2))
                min_nz = int(pm.group(3))
            break

    if start_index is None:
        raise ValueError("Could not find RESNUCLEI matrix header.")
    if None in (max_z, max_nz, min_nz, k):
        raise ValueError("Could not parse matrix parameters (max_z, max_nz, min_nz, k).")

    # Find end of numeric block (next starred separator)
    for i, line in enumerate(lines[start_index:]):
        if line.strip().startswith("*"):
            end_index = start_index + i
            break
    if end_index is None:
        end_index = len(lines)

    # Parse all scientific-notation numbers from the block
    block = "\n".join(lines[start_index:end_index])
    values = [float(x) for x in re.findall(r"[-+]?\d+\.\d+E[-+]\d+", block)]

    m_count = max_nz - min_nz + 1
    expected = max_z * m_count
    if len(values) < expected:
        raise ValueError(f"Not enough matrix values: got {len(values)}, expected {expected}")

    # IMPORTANT: Fortran ordering (Z index varies fastest in stream)
    raw = np.array(values[:expected], dtype=float).reshape((max_z, m_count), order="F")

    arr = np.asarray(raw, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"rndata must be 2D, got shape {arr.shape}")

    z_count, m_count = arr.shape
    z_values = np.arange(1, 1 + z_count, dtype=int)
    m_values = np.arange(1, 1 + m_count, dtype=int)

    n_min = np.max([int(z_values[0] + m_values[0] + k),0])
    n_max = int(z_values[-1] + m_values[-1] + k)
    n_values = np.arange(n_min, n_max + 1, dtype=int)

    out = np.zeros((z_count, n_values.size), dtype=float)

    for iz, z in enumerate(z_values):
        for jm, m in enumerate(m_values):
            n = int(z + m + k)
            out[iz, n - n_min - 1] = arr[iz, jm]

    return RESNCULEi_output(z_axis=z_values, n_axis=n_values, values=out)

def combine_outputs(outputs: list[RESNCULEi_output]):
    # Get global z and n axes
    all_z = np.unique(np.concatenate([o.z_axis for o in outputs]))
    all_n = np.unique(np.concatenate([o.n_axis for o in outputs]))

    combined_values = np.zeros((all_z.size, all_n.size), dtype=float)

    for o in outputs:
        z_indices = np.searchsorted(all_z, o.z_axis)
        n_indices = np.searchsorted(all_n, o.n_axis)

        for iz, z_idx in enumerate(z_indices):
            for in_, n_idx in enumerate(n_indices):
                combined_values[z_idx, n_idx] += o.values[iz, in_]
    combined_values /= len(outputs)  # Average if needed

    return RESNCULEi_output(z_axis=all_z, n_axis=all_n, values=combined_values)

file_in_folder = RESNUCLEI_OUTPUT_FILE

def calculate_isotope_production_per_folder(folder):
    file = Path(folder) / file_in_folder
    if not file.exists():
        raise FileNotFoundError(f"Expected RESNUCLEi output file not found: {file}")
    return parse_resnuclei_output(file)

def calculate_isotope_production_per_project(folders):
    outputs = []
    for folder in folders:
        try:
            output = calculate_isotope_production_per_folder(folder)
            outputs.append(output)
        except Exception as e:
            print(f"Error processing {folder}: {e}")
    if not outputs:
        raise ValueError("No valid RESNUCLEi outputs found for project.")
    combined = combine_outputs(outputs)
    return combined

def process_isotope_production_material(project_folders, output_file=None):
    project_isotopes = calculate_isotope_production_per_project(project_folders)
    if output_file:
        with open(output_file, "w") as f:
            json.dump({"z_axis": project_isotopes.z_axis.tolist(), "n_axis": project_isotopes.n_axis.tolist(), "values": project_isotopes.values.tolist()}, f)
    return project_isotopes

############## Main Processing Functions #############

def post_proc(project_folders):
    output_folder = DEFAULT_OUTPUT_ROOT
    output = {}
    for material, folders in project_folders.items():
        yields = process_neutron_yield_material(folders, output_file=output_folder / f"{material}_neutron_yield.json")
        isotopes = process_isotope_production_material(folders, output_file=output_folder / f"{material}_isotope_production.json")
        output[material] = {
            "yields": yields,
            "isotopes": isotopes
        }
    return output