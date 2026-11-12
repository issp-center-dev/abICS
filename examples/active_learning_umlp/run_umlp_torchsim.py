#!/usr/bin/env python3
"""Relax abICS reference structures in GPU batches with TorchSim and Orb.

Pinned dependencies and the complete workflow are documented in README.md.
Run once per reference batch, in place of parallel_run.sh.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import os
from pathlib import Path
import time

import ase.io
import numpy as np
import torch
import torch_sim as ts
from orb_models.forcefield import pretrained
from torch_sim.models.orb import OrbModel
from torch_sim.optimizers.cell_filters import CellFilter


class MeasuredOrbModel(OrbModel):
    """Record actual model batch sizes, including batches shrunk by convergence."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.max_batch_seen = 0
        self.forward_calls = 0

    def forward(self, state, **kwargs):
        self.max_batch_seen = max(self.max_batch_seen, state.n_systems)
        self.forward_calls += 1
        return super().forward(state, **kwargs)


def load_model(device="cuda", precision="float32-highest"):
    dtype = torch.float64 if precision == "float64" else torch.float32
    options = {}
    if local_weights := os.environ.get("ORB_WEIGHTS_PATH"):
        options["weights_path"] = str(Path(local_weights).expanduser())
    forcefield, adapter = pretrained.orb_v3_conservative_20_omat(
        device=device, precision=precision, compile=False, **options,
    )
    return MeasuredOrbModel(forcefield, adapter, device=device, dtype=dtype)


def synchronize(model):
    if model.device.type == "cuda":
        torch.cuda.synchronize(model.device)


def relax_batch(atoms, model, *, fire_steps=200, max_steps=1000, fmax=0.04):
    """Return states in input order; optimize positions AND the full cell."""
    system = ts.io.atoms_to_state(atoms, device=model.device, dtype=model.dtype)
    if fire_steps:
        system = ts.optimize(
            system, model, optimizer=ts.Optimizer.fire,
            convergence_fn=ts.generate_force_convergence_fn(
                force_tol=0.2, include_cell_forces=True,
            ),
            max_steps=fire_steps, steps_between_swaps=1, autobatcher=False,
            init_kwargs={"cell_filter": CellFilter.frechet},
        )
        # A FIRE OptimState must not be passed directly to the BFGS runner:
        # the runner would skip BFGS initialization for an existing OptimState.
        system = ts.SimState.from_state(system)
    convergence = ts.generate_force_convergence_fn(
        force_tol=fmax, include_cell_forces=True,
    )
    final = ts.optimize(
        system, model, optimizer=ts.Optimizer.bfgs,
        convergence_fn=convergence, max_steps=max_steps,
        steps_between_swaps=1, autobatcher=False,
        init_kwargs={"cell_filter": CellFilter.frechet},
    )
    # TorchSim returns states even when max_steps is reached: check explicitly.
    converged = convergence(final).detach().cpu().tolist()
    atomic_fmax = ts.system_wise_max_force(final).detach().cpu().tolist()
    cell_fmax = final.cell_forces.norm(dim=2).amax(dim=1).detach().cpu().tolist()
    # Reevaluate final energies at the returned positions/cells.
    prediction = model(ts.SimState.from_state(final))
    energies = prediction["energy"].detach().cpu().tolist()
    structures = ts.io.state_to_atoms(final)
    return structures, [
        {"converged": bool(ok), "energy_eV": float(energy),
         "atomic_fmax": float(af), "cell_fmax": float(cf)}
        for ok, energy, af, cf in zip(converged, energies, atomic_fmax, cell_fmax, strict=True)
    ]


def write_json(path, value):
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    os.replace(temp, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rundirs", type=Path, default=Path("rundirs.txt"))
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--precision", default="float32-highest",
                        choices=["float32-high", "float32-highest", "float64"])
    parser.add_argument("--fire-steps", type=int, default=200)
    parser.add_argument("--max-steps", type=int, default=1000)
    parser.add_argument("--fmax", type=float, default=0.04)
    parser.add_argument("--summary", type=Path, default=Path("torchsim-summary.json"))
    args = parser.parse_args()
    if args.batch_size < 1 or args.fire_steps < 0 or args.max_steps < 1 or not math.isfinite(args.fmax) or args.fmax <= 0:
        parser.error("require batch-size/max-steps > 0, fire-steps >= 0, fmax > 0")
    dirs = [Path(line.strip()).resolve() for line in args.rundirs.read_text().splitlines()
            if line.strip()]
    if not dirs or len(set(dirs)) != len(dirs):
        parser.error("rundirs must contain unique calculation directories")
    # Validate the entire list before starting; avoid mixing results from runs.
    for directory in dirs:
        if not (directory / "structure.vasp").is_file():
            parser.error(f"missing structure.vasp: {directory}")
        if any((directory / name).exists() for name in
               ("energy.dat", "structure_norel.vasp", "torchsim_relaxation.json")):
            parser.error(f"existing reference outputs: {directory}; use fresh input directories")
    model = load_model(args.device, args.precision)
    if model.device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(model.device)
    summary = {
        "model": "orb_v3_conservative_20_omat", "precision": args.precision,
        "cell_filter": "frechet", "batch_size": args.batch_size,
        "fire_steps": args.fire_steps, "max_steps": args.max_steps, "fmax": args.fmax,
        "versions": {p: importlib.metadata.version(p) for p in
                     ("torch", "torch-sim-atomistic", "orb-models", "ase")},
        "batches": [], "structures": [],
    }
    for offset in range(0, len(dirs), args.batch_size):
        chunk = dirs[offset:offset + args.batch_size]
        atoms = [ase.io.read(d / "structure.vasp") for d in chunk]
        synchronize(model)
        start = time.perf_counter()
        relaxed, records = relax_batch(atoms, model, fire_steps=args.fire_steps,
                                      max_steps=args.max_steps, fmax=args.fmax)
        synchronize(model)
        elapsed = time.perf_counter() - start
        for directory, original, final, record in zip(chunk, atoms, relaxed, records, strict=True):
            if not np.array_equal(original.numbers, final.numbers):
                raise RuntimeError(f"atomic ordering changed: {directory}")
            finite = (math.isfinite(record["energy_eV"])
                      and np.isfinite(final.positions).all()
                      and np.isfinite(final.cell.array).all()
                      and final.get_volume() > 0
                      and math.isfinite(record["atomic_fmax"])
                      and math.isfinite(record["cell_fmax"]))
            if not finite:
                raise RuntimeError(f"Non-finite relaxation output or invalid cell: {directory}")
            record.update(directory=str(directory), converged=bool(record["converged"]))
            ase.io.write(directory / "structure_norel.vasp", original, format="vasp")
            if record["converged"]:
                temp = directory / "structure.vasp.tmp"
                ase.io.write(temp, final, format="vasp")
                os.replace(temp, directory / "structure.vasp")
                temp = directory / "energy.dat.tmp"
                temp.write_text(f'{record["energy_eV"]:.16g}\n')
                os.replace(temp, directory / "energy.dat")
            else:
                ase.io.write(directory / "structure_unconverged.vasp", final, format="vasp")
            write_json(directory / "torchsim_relaxation.json", record)
            summary["structures"].append(record)
        summary["batches"].append({"n_structures": len(chunk), "relaxation_seconds": elapsed})
        summary.update(max_batch_seen=model.max_batch_seen, forward_calls=model.forward_calls,
                       peak_gpu_allocated_bytes=(torch.cuda.max_memory_allocated(model.device)
                                                 if model.device.type == "cuda" else 0))
        write_json(args.summary, summary)
        print(f"Finished {offset + len(chunk)}/{len(dirs)} structures in {elapsed:.3f} s; "
              f"converged {sum(r['converged'] for r in records)}/{len(records)}", flush=True)
    failed = [r["directory"] for r in summary["structures"] if not r["converged"]]
    if failed:
        raise SystemExit(f"{len(failed)} unconverged references; no energy.dat written for them. "
                         f"Inspect {args.summary} before collecting results.")


if __name__ == "__main__":
    main()
