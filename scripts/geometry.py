from __future__ import annotations
from pathlib import Path

import numpy as np
import pyg4ometry as pg4
from pygeomhpges.materials import make_enriched_germanium
from pygeomtools.materials import LegendMaterialRegistry
import tempfile
import os


####################### SIMULATION #######################

HEIGHT_IN_G_OVER_CM2 = 2000
RADIUS_IN_CM = 125
DENSITIES = {
    "lar": 1.396,  # g/cm^3
    "water": 1.0,  # g/cm^3
    "rock": 2.65,  # g/cm^3
}
AVERAGE_A = {
    "lar": 40.0,
    "water": 14.335,
    "rock": 22.0,
}

def geometry(mat_name: str):
    reg = pg4.geant4.Registry()
    matreg = LegendMaterialRegistry(reg, enable_optical=False)

    if mat_name == "lar":
        mat = matreg.liquidargon
    elif mat_name == "water":
        mat = matreg.water
    elif mat_name == "rock":
        mat = matreg.rock
    elif mat_name == "enrGe":
        mat = make_enriched_germanium(reg, "enrGe", enrichment=0.92)
    else:
        msg = f"unknown material {mat_name}"
        raise ValueError(msg)

    height = HEIGHT_IN_G_OVER_CM2 / DENSITIES[mat_name]  # convert to cm

    world_s = pg4.geant4.solid.Tubs(
        "world",
        0,
        RADIUS_IN_CM,
        height,
        0,
        2 * np.pi,
        registry=reg,
        lunit="cm",
    )
    world_l = pg4.geant4.LogicalVolume(world_s, mat, "world", registry=reg)
    reg.setWorld(world_l)

    safety = 1e-5
    detector_s = pg4.geant4.solid.Tubs(
        "detector",
        0,
        RADIUS_IN_CM - safety,
        height - safety,
        0,
        2 * np.pi,
        registry=reg,
        lunit="cm",
    )
    detector_l = pg4.geant4.LogicalVolume(detector_s, mat, "detector", registry=reg)
    pg4.geant4.PhysicalVolume(
        [0, 0, 0], [0, 0, 0], detector_l, "detector", world_l, registry=reg
    )

    return reg

def patch_fluka_body_store_keyerror():
    """Patch older pyg4ometry body-store implementations that can raise KeyError.

    Some pyg4ometry builds use `return self.hashBody[body.hash()]` directly.
    If internal maps get out-of-sync, conversion crashes. This patch uses the
    safer behavior: return cached body if available, otherwise add and return.
    """
    try:
        from pyg4ometry.fluka import fluka_registry as fr
    except Exception:
        return 0

    def _safe_get_degenerate_body(self, body):
        body_hash = body.hash()
        hash_body = getattr(self, "hashBody", None)
        if isinstance(hash_body, dict) and body_hash in hash_body:
            return hash_body[body_hash]

        add_body = getattr(self, "addBody", None)
        if callable(add_body):
            return add_body(body)

        # Conservative fallback if this store has no addBody method.
        if isinstance(hash_body, dict):
            hash_body[body_hash] = body
        name_body = getattr(self, "nameBody", None)
        if isinstance(name_body, dict) and hasattr(body, "name"):
            name_body[body.name] = body
        return body

    patched = 0
    for cls_name in ("FlukaBodyStore", "FlukaBodyStoreExact"):
        cls = getattr(fr, cls_name, None)
        if cls is None:
            continue
        if hasattr(cls, "getDegenerateBody"):
            cls.getDegenerateBody = _safe_get_degenerate_body
            patched += 1

    return patched


def patch_elliptical_tube_radius_bug():
    """Patch pyg4ometry elliptical tube conversion for correct REC semi-axes.

    GDML/Geant4 `EllipticalTube` uses semi-axes (half-lengths) for dx/dy/dz.
    Some pyg4ometry versions divide dx/dy by 2 again when building FLUKA REC,
    which halves the tube radii in the output.
    """
    try:
        from pyg4ometry.convert import geant42Fluka as g2f
    except Exception:
        return False

    def _fixed_geant4_elliptical_tube_2_fluka(
        flukaNameCount,
        solid,
        mtra=g2f._np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        tra=g2f._np.array([0, 0, 0]),
        flukaRegistry=None,
        addRegistry=True,
        commentName="",
        bakeTransform=False,
    ):
        name = format(flukaNameCount, "04")

        from pyg4ometry.gdml import Units as _Units

        rotation = g2f._transformation.matrix2tbxyz(mtra)
        transform = g2f._rotoTranslationFromTra2(
            "T" + name, [rotation, tra], flukaregistry=flukaRegistry
        )

        uval = _Units.unit(solid.lunit) / 10
        pDx = solid.evaluateParameter(solid.pDx) * uval
        pDy = solid.evaluateParameter(solid.pDy) * uval
        pDz = solid.evaluateParameter(solid.pDz) * uval

        # pDx/pDy/pDz are Geant4 half-lengths (semi-axes / half-height).
        # FLUKA REC needs semi-axes for the ellipse and full axis vector length.
        if not bakeTransform:
            fbody1 = g2f._fluka.REC(
                "B" + name + "01",
                [0, 0, -pDz],
                [0, 0, 2 * pDz],
                [pDx, 0, 0],
                [0, pDy, 0],
                transform=transform,
                flukaregistry=flukaRegistry,
                comment=commentName,
            )
        else:
            fbody1 = g2f._fluka.REC(
                "B" + name + "01",
                mtra @ g2f._np.array([0, 0, -pDz]) + tra / 10,
                mtra @ g2f._np.array([0, 0, 2 * pDz]),
                mtra @ g2f._np.array([pDx, 0, 0]),
                mtra @ g2f._np.array([0, pDy, 0]),
                transform=None,
                flukaregistry=flukaRegistry,
                comment=commentName,
            )

        fzone = g2f._fluka.Zone()
        fzone.addIntersection(fbody1)

        fregion = g2f._fluka.Region("R" + name)
        fregion.addZone(fzone)

        flukaNameCount += 1
        return fregion, flukaNameCount

    g2f.geant4EllipticalTube2Fluka = _fixed_geant4_elliptical_tube_2_fluka
    return True


def sanitize_material_cards(path):
    """Fix known pyg4ometry MATERIAL card field placement issues per FLUKA manual.

    MATERIAL card format: Z  [atomic_weight_override]  density  [mat_num]  [alt_mat]  [mass_number]  SDUM
    
    - `* element-simple` cards: WHAT(2) blank (FLUKA infers), WHAT(6)=0 (natural isotopic composition)
    - `* isotope` cards: WHAT(6) = isotope mass number (integer)
    """
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()

    changed = 0
    context = None

    for i, line in enumerate(lines):
        stripped = line.strip()

        if stripped.startswith("* element-simple:"):
            context = "element-simple"
            continue

        if stripped.startswith("* isotope:"):
            context = "isotope"
            continue

        if not stripped.startswith("MATERIAL"):
            continue

        parts = [p.strip() for p in line.rstrip("\n").split(",")]
        if len(parts) < 8:
            parts.extend([""] * (8 - len(parts)))

        # element-simple: WHAT(2) blank, WHAT(6)=0 for natural composition
        if context == "element-simple":
            if parts[2] != "":
                parts[2] = ""
                changed += 1
            if parts[6] != "0" and parts[6] != "":
                parts[6] = "0"
                changed += 1

        # isotope: WHAT(2) blank, WHAT(6) = mass number (integer)
        # If WHAT(6) is non-integer, clear it (isotope handling may need COMPOUND card)
        if context == "isotope":
            if parts[2] != "":
                parts[2] = ""
                changed += 1

        lines[i] = ", ".join(parts) + "\n"
        context = None

    if changed:
        with open(path, "w", encoding="utf-8") as fh:
            fh.writelines(lines)

    return changed


def convert_gdml_to_fluka(gdml, out):
    reader = pg4.gdml.Reader(gdml)
    greg = reader.getRegistry()

    # First pass: default conversion path.
    try:
        freg = pg4.convert.geant4Reg2FlukaReg(greg)
    except KeyError as exc:
        print(f"Primary conversion failed with KeyError: {exc}")
        print("Retrying with bakeTransforms=True")
        freg = pg4.convert.geant4Reg2FlukaReg(greg, bakeTransforms=True)

    w = pg4.fluka.Writer()
    w.addDetector(freg)
    w.write(out)

def strip_leading_free_card(geometry_text: str) -> str:
	lines = geometry_text.splitlines()
	if lines and lines[0].strip().upper() == "FREE":
		lines = lines[1:]
	return "\n".join(lines).strip()


def generate_geometry(mat_name: str):
    with tempfile.TemporaryDirectory() as tmpdirname:
        gdml_path = os.path.join(tmpdirname, "geometry.gdml")
        fluka_path = os.path.join(tmpdirname, "geometry.inp")

        reg = geometry(mat_name)
        z_extent = reg.logicalVolumeDict["world"].solid.pDz
        w = pg4.gdml.Writer()
        w.addDetector(reg)
        w.write(gdml_path)

        convert_gdml_to_fluka(gdml_path, fluka_path)
        sanitize_material_cards(fluka_path)

        text = Path(fluka_path).read_text(encoding="utf-8")

        return {
            "text": strip_leading_free_card(text),
            "z_extent_cm": z_extent
        }