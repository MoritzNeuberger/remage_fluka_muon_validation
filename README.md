# Muon simulations in FLUKA for validation of remage

This repository contains a small workflow to run muon simulations in [FLUKA](https://fluka.cern/) for simple geometry with `lar`, `water`, and `rock` materials. The purpose of this is to generate reference observables for validating the [remage](https://github.com/legend-exp/remage) simulation framework based on [Geant4](https://geant4.web.cern.ch/). A comparison is shown in the validation suite.

The workflow will generate FLUKA input files, compile and run them, and post-process neutron yield and isotope production outputs into JSON files. To run the full workflow, the user has to provide FLUKA binaries and set up a Python environment.

## 1) Setup

### Python environment

From this folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### FLUKA binaries location

For default behavior, FLUKA binaries must exist at:

```text
{repository root}/sw/fluka4-5.1/bin/
```

and this directory must contain at least:

- `fff`
- `lfluka`
- `rfluka`

Why: project generation writes a per-project `compile_and_run.sh` that calls those tools.

The default path is set in `scripts/project.py` via `--fluka-bin-path` default:

```text
sw/fluka4-5.1/bin
```


## 2) Default use

Run the full default workflow (generate projects, run FLUKA, post-process):

```bash
cd scripts
python3 orchestrate.py
```

Default values are:

- 20 projects per material
- 20 threads
- 100000 primaries per project
- materials: `lar`, `water`, `rock`

Output JSON files are written to:

```text
gen/{material}_neutron_yield.json
gen/{material}_isotope_production.json
```

You can reduce runtime for a quick smoke test:

```bash
python3 orchestrate.py --n-projects-per-material 1 --n-threads 1 --n-primaries 1000
```


## 3) Where geometry is defined

Geometry for this workflow is defined in:

- `scripts/geometry.py`

Main places to edit:

- `HEIGHT_IN_G_OVER_CM2`
- `RADIUS_IN_CM`
- `DENSITIES`
- `geometry(mat_name)` (the actual Geant4 solids and materials)

Current implementation uses cylindrical `Tubs` for world and detector volumes.


## 4) How to align geometry with the latest remage validation

remage validation reference link: TODO (insert link here).

1. Update material definitions and densities in `scripts/geometry.py` to match the validation reference.
2. Update dimensions in `HEIGHT_IN_G_OVER_CM2`, `RADIUS_IN_CM`, and solid construction in `geometry(mat_name)`.
3. Keep units consistent:
	- density in g/cm^3
	- geometry lengths in cm
4. Regenerate one project and verify `fluka_input.inp` geometry/cards are coherent.

Notes:

- The z-extent from geometry is used downstream in `scripts/project.py` to set `BEAMPOS` and `USRICALL` values.
- If geometry changes, those beam/cut values update automatically because they are derived from `z_extent_cm`.


## 5) Useful direct commands

Generate one project folder only (no run):

```bash
cd scripts
python3 project.py --material lar
```

Run a full orchestrated job with custom scale:

```bash
cd scripts
python3 orchestrate.py --n-projects-per-material 5 --n-threads 5 --n-primaries 50000
```

