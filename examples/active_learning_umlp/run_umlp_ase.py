#!/usr/bin/env python3
"""Relax one reference structure with ASE and Orb 0.7 (CPU or CUDA)."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import ase.io
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS, FIRE
from orb_models.forcefield import pretrained
from orb_models.forcefield.inference.calculator import ORBCalculator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    directory = args.directory.resolve()
    if not (directory / "structure.vasp").is_file():
        parser.error(f"missing structure.vasp: {directory}")
    if any((directory / name).exists() for name in
           ("energy.dat", "structure_norel.vasp", "ase_relaxation.json")):
        parser.error("existing reference outputs; use fresh input directories")
    options = {}
    if local_weights := os.environ.get("ORB_WEIGHTS_PATH"):
        options["weights_path"] = str(Path(local_weights).expanduser())
    model, adapter = pretrained.orb_v3_conservative_20_omat(
        device=args.device, precision="float32-highest", compile=False, **options,
    )
    atoms = ase.io.read(directory / "structure.vasp")
    ase.io.write(directory / "structure_norel.vasp", atoms, format="vasp")
    atoms.calc = ORBCalculator(model, adapter, device=args.device)
    filtered = FrechetCellFilter(atoms)
    FIRE(filtered, logfile=str(directory / "fire.log")).run(fmax=0.2, steps=200)
    converged = BFGS(filtered, logfile=str(directory / "bfgs.log")).run(
        fmax=0.04, steps=1000,
    )
    energy = float(atoms.get_potential_energy())
    result = {"converged": bool(converged and math.isfinite(energy)),
              "energy_eV": energy}
    (directory / "ase_relaxation.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    if not result["converged"]:
        ase.io.write(directory / "structure_unconverged.vasp", atoms, format="vasp")
        raise SystemExit("BFGS did not converge; no energy.dat written")
    temporary = directory / "structure.vasp.tmp"
    ase.io.write(temporary, atoms, format="vasp")
    temporary.replace(directory / "structure.vasp")
    temporary = directory / "energy.dat.tmp"
    temporary.write_text(f"{energy:.16g}\n")
    temporary.replace(directory / "energy.dat")


if __name__ == "__main__":
    main()
