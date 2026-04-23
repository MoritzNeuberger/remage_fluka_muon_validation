import json
import awkward as ak
import numpy as np
from dataclasses import dataclass
import re
from pygeomhpges.materials import enriched_germanium_density
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[1]
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "gen"
NEUTRON_SCORING_FILE = "fluka_input001_neutron_scoring.dat"
FLUKA_INPUT_FILE = "fluka_input.inp"
RESNUCLEI_OUTPUT_FILE = "fluka_input001.out"


def _mode_folder(folder, mode):
    """Return mode subfolder when available, otherwise keep legacy folder layout."""
    folder = Path(folder)
    candidate = folder / mode
    if candidate.exists() and candidate.is_dir():
        return candidate
    return folder

DENSITIES = {
    "lar": 1.396,  # g/cm^3
    "water": 1.0,  # g/cm^3
    "rock": 2.65,  # g/cm^3
    "enrGe": enriched_germanium_density(0.92).magnitude,  # g/cm^3 (approximate for enriched germanium)
}
ELEMENTS_PER_MATERIAL = {
    "lar": [("Ar", 18)],
    "water": [("H", 1), ("O", 8)],
    "rock": [("O", 8), ("Si", 14)],  # O, Si
    "enrGe": [("Ge", 32)],
}

def get_n_muons_per_folder(folder):
    folder = Path(folder)
    with (folder / NEUTRON_SCORING_FILE).open() as f:
        lines = f.readlines()
    data = []
    for line in lines:        
        if line.startswith("#"):
            continue
        line_split = line.strip().split()
        try:
            data.append(int(line_split[0]))
        except (ValueError, IndexError):
            continue
    return len(data)

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
    material = str(Path(folder)).split("/")[-2].split("_")[0].split("-")[0]
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
    yields = []
    for folder in folders:
        neutron_folder = _mode_folder(folder, "neutron")
        #try:
        yields.append(calculate_neutron_yield_per_folder(neutron_folder))
        #except Exception as exc:
        #    print(f"[warn] skipping neutron yield for {neutron_folder}: {exc}")

    if not yields:
        raise ValueError("No valid neutron yield outputs found for project.")

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
    z_values = np.arange(1, z_count + 1, dtype=int)
    m_values = np.arange(1, m_count + 1, dtype=int)

    n_min = max(int(z_values[0] + m_values[0] + k), 0)
    n_max = int(z_values[-1] + m_values[-1] + k)
    n_values = np.arange(n_min, n_max + 1, dtype=int)

    out = np.zeros((z_count, n_values.size), dtype=float)

    for iz, z in enumerate(z_values):
        for jm, m in enumerate(m_values):
            n = int(z + m + k)
            if n < n_min or n > n_max:
                continue
            out[iz, n - n_min] = arr[iz, jm]

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
            isotope_folder = _mode_folder(folder, "isotope")
            output = calculate_isotope_production_per_folder(isotope_folder)
            outputs.append(output)
        except Exception as e:
            print(f"[warn] skipping isotope output for {folder}: {e}")
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

############## Neutron Production Process Analysis #############
fluka_pid_lookup = {
    -6: {"pdg":	9999, "name":	"4-HELIUM"},
    -5: {"pdg":	9999, "name":	"3-HELIUM"},
    -4: {"pdg":	9999, "name":	"TRITON"},
    -3: {"pdg":	9999, "name":	"DEUTERON"},
    -2: {"pdg":	9999, "name":	"HEAVYION"},
    -1: {"pdg":	9999, "name":	"OPTIPHOT"},
    0: {"pdg":	9999, "name":	"RAY"},
    1: {"pdg":	2212, "name":	"PROTON"},
    2: {"pdg":	-2212, "name":	"APROTON"},
    3: {"pdg":	11, "name":	"ELECTRON"},
    4: {"pdg":	-11, "name":	"POSITRON"},
    5: {"pdg":	12, "name":	"NEUTRIE"},
    6: {"pdg":	-12, "name":	"ANEUTRIE"},
    7: {"pdg":	22, "name":	"PHOTON"},
    8: {"pdg":	2112, "name":	"NEUTRON"},
    9: {"pdg":	-2112, "name":	"ANEUTRON"},
    10: {"pdg":	-13, "name":	"MUON+"},
    11: {"pdg":	13, "name":	"MUON-"},
    12: {"pdg":	130, "name":	"KAONLONG"},
    13: {"pdg":	211, "name":	"PION+"},
    14: {"pdg":	-211, "name":	"PION-"},
    15: {"pdg":	321, "name":	"KAON+"},
    16: {"pdg":	-321, "name":	"KAON-"},
    17: {"pdg":	3122, "name":	"LAMBDA"},
    18: {"pdg":	-3122, "name":	"ALAMBDA"},
    19: {"pdg":	310, "name":	"KAONSHRT"},
    20: {"pdg":	3112, "name":	"SIGMA-"},
    21: {"pdg":	3222, "name":	"SIGMA+"},
    22: {"pdg":	3212, "name":	"SIGMAZER"},
    23: {"pdg":	111, "name":	"PIZERO"},
    24: {"pdg":	311, "name":	"KAONZERO"},
    25: {"pdg":	-311, "name":	"AKAONZER"},
    26: {"pdg":	0, "name":	"RESERVED"},
    27: {"pdg":	14, "name":	"NEUTRIM"},
    28: {"pdg":	-14, "name":	"ANEUTRIM"},
    29: {"pdg":	0, "name":	"RESERVED"},
    30: {"pdg":	0, "name":	"RESERVED"},
    31: {"pdg":	-3222, "name":	"ASIGMA-"},
    32: {"pdg":	-3212, "name":	"ASIGMAZE"},
    33: {"pdg":	-3112, "name":	"ASIGMA+"},
    34: {"pdg":	3322, "name":	"XSIZERO"},
    35: {"pdg":	-3322, "name":	"AXSIZERO"},
    36: {"pdg":	3312, "name":	"XSI-"},
    37: {"pdg":	-3312, "name":	"AXSI+"},
    38: {"pdg":	3334, "name":	"OMEGA-"},
    39: {"pdg":	-3334, "name":	"AOMEGA+"},
    40: {"pdg":	2112, "name":	"WWLOWNEU"},
    41: {"pdg":	-15, "name":	"TAU+"},
    42: {"pdg":	15, "name":	"TAU-"},
    43: {"pdg":	16, "name":	"NEUTRIT"},
    44: {"pdg":	-16, "name":	"ANEUTRIT"},
    45: {"pdg":	411, "name":	"D+"},
    46: {"pdg":	-411, "name":	"D-"},
    47: {"pdg":	421, "name":	"D0"},
    48: {"pdg":	-421, "name":	"D0BAR"},
    49: {"pdg":	431, "name":	"DS+"},
    50: {"pdg":	-431, "name":	"DS-"},
    51: {"pdg":	4122, "name":	"LAMBDAC+"},
    52: {"pdg":	4232, "name":	"XSIC+"},
    53: {"pdg":	4132, "name":	"XSIC0"},
    54: {"pdg":	4322, "name":	"XSIPC+"},
    55: {"pdg":	4312, "name":	"XSIPC0"},
    56: {"pdg":	4332, "name":	"OMEGAC0"},
    57: {"pdg":	-4122, "name":	"ALAMBDC-"},
    58: {"pdg":	-4232, "name":	"AXSIC-"},
    59: {"pdg":	-4132, "name":	"AXSIC0"},
    60: {"pdg":	-4322, "name":	"AXSIPC-"},
    61: {"pdg":	-4312, "name":	"AXSIPC0"},
    62: {"pdg":	-4332, "name":	"AOMEGAC0"},
    63: {"pdg":	9999, "name":	"RESERVED"},
    64: {"pdg":	9999, "name":	"RESERVED"}
}


def file_parser(lines):
    data = []
    process_current = {"primary": None, "secundaries": {}, "n_secundaries": None, "proc_info": {}}
    for line in lines:
        if line.startswith("#") or line.startswith("*"):
            continue
        if line.startswith("I"):
            #'I            2  101      1    7      3  1.74676E-01 -4.97879E-01 -5.50250E+02 F\n',
            if process_current["primary"] is not None:
                process_current["proc_info"]["is_primary_in_secundaries"] = process_current["primary"] in process_current["secundaries"].keys()
                data.append(process_current.copy())
            process_current =  {"primary": None, "secundaries": {}, "n_secundaries": None, "proc_info": {}}
            
            if "\\x00" in line:
                line = line.replace("\\x00", "")
            lines_split = line.split()
            try:
                process_current["primary"] = fluka_pid_lookup[int(lines_split[4])]["pdg"]
            except (KeyError, ValueError):
                continue
            process_current["proc_info"]["proc_id"] = int(lines_split[2])
            process_current["n_secundaries"] = int(lines_split[5])

        if line.startswith("O"):
            #'O            2      1    8 T\n',
            lines_split = line.split()
            if "*" in lines_split[3]:
                continue
            if fluka_pid_lookup[int(lines_split[3])]["pdg"] not in process_current["secundaries"].keys():
                process_current["secundaries"][fluka_pid_lookup[int(lines_split[3])]["pdg"]] = 0
            process_current["secundaries"][fluka_pid_lookup[int(lines_split[3])]["pdg"]] += 1
    return data

def calculate_neutron_multiplicity_per_primary(proc_data, primary):
    multiplicities = []
    for process in proc_data:
        if process["primary"] == primary:
            mult = process["secundaries"].get(2112, 0)
            if primary == 2112 and mult > 0:
                mult -= 1
            multiplicities.append(mult)
    return multiplicities

def calculate_neutron_multiplicity_per_folder(folder):

    file = Path(folder) / "fluka_input001_neutron_interaction_particles.dat"

    with Path(file).open() as f:
        lines = f.readlines()
    proc_data = file_parser(lines)

    neutron_mults_per_primary = {}
    unique_primaries = set([process["primary"] for process in proc_data])
    for prim in unique_primaries:
        neutron_mults_per_primary[prim] = calculate_neutron_multiplicity_per_primary(proc_data, prim)
    for prim in neutron_mults_per_primary.keys():
        neutron_mults_per_primary[prim] = np.bincount(neutron_mults_per_primary[prim], minlength=25)

    n_muon = get_n_muons_per_folder(folder)
    for prim in neutron_mults_per_primary.keys():
        neutron_mults_per_primary[prim] = neutron_mults_per_primary[prim] / n_muon
    return neutron_mults_per_primary

def calculate_neutron_multiplicity_per_project(folders):
    all_primaries = set()
    folder_mults_list = []
    for folder in folders:
        neutron_folder = _mode_folder(folder, "neutron")
        try:
            folder_mults = calculate_neutron_multiplicity_per_folder(neutron_folder)
            folder_mults_list.append(folder_mults)
            all_primaries.update(folder_mults.keys())
        except Exception as exc:
            print(f"[warn] skipping neutron multiplicity for {neutron_folder}: {exc}")

    if not folder_mults_list:
        raise ValueError("No valid neutron multiplicity outputs found for project.")

    combined_mults = {prim: [] for prim in all_primaries}
    for folder_mults in folder_mults_list:
        for prim in all_primaries:
            if prim in folder_mults:
                combined_mults[prim].append(folder_mults[prim])
            else:
                combined_mults[prim].append(np.zeros(25))
    averaged_mults = {prim: np.mean(combined_mults[prim], axis=0).tolist() for prim in all_primaries}
    return averaged_mults


def process_neutron_multiplicity(project_folders, output_file=None):
    project_neutron_mults = calculate_neutron_multiplicity_per_project(project_folders)
    if output_file:
        with open(output_file, "w") as f:
            json.dump(project_neutron_mults, f)
    return project_neutron_mults

############## Main Processing Functions #############

def post_proc(project_folders):
    output_folder = DEFAULT_OUTPUT_ROOT
    output = {}
    for case, folders in project_folders.items():
        try:
            yields = process_neutron_yield_material(
                folders,
                output_file=output_folder / f"{case[0]}-{case[1]}GeV_neutron-yield.json",
            )
        except Exception as exc:
            print(f"[warn] neutron yield aggregation failed for {case}: {exc}")
            yields = None

        try:
            isotopes = process_isotope_production_material(
                folders,
                output_file=output_folder / f"{case[0]}-{case[1]}GeV_isotope-production.json",
            )
        except Exception as exc:
            print(f"[warn] isotope aggregation failed for {case}: {exc}")
            isotopes = None

        try:
            neutron_mult = process_neutron_multiplicity(
                folders,
                output_file=output_folder / f"{case[0]}-{case[1]}GeV_neutron-multiplicity.json",
            )
        except Exception as exc:
            print(f"[warn] neutron multiplicity aggregation failed for {case}: {exc}")
            neutron_mult = None

        output[case] = {
            "yields": yields,
            "isotopes": isotopes,
            "neutron_mult": neutron_mult
        }
    return output