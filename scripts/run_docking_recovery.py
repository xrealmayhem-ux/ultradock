#!/usr/bin/env python3
"""Prepare localized 5-HT2A receptor models and run a bounded validation pilot."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "docking_recovery"
AUDIT = ROOT / "results" / "docking_audit"
VINA_GPU_DIR = Path(
    os.environ.get("ULTRADOCK_VINA_GPU_DIR", Path.home() / "Vina-GPU")
).expanduser().resolve()
VINA_GPU_EXECUTABLE = VINA_GPU_DIR / "Vina-GPU"
VINA_GPU_WRAPPER = VINA_GPU_DIR / "run_vina.sh"
VINA_GPU_KERNEL = VINA_GPU_DIR / "Kernel2_Opt.bin"
VINA_GPU_PARAMETERS = {
    "scoring": "vina",
    "thread": 1000,
    "search_depth": 0,
    "num_modes": 9,
    "energy_range": 3,
}
LIGAND_SDF = ROOT / "structures" / "ligands" / "sdf" / "docking_ligands_selected_3d.sdf"
SELECTION = ROOT / "data" / "processed" / "ligand_selection.csv"
ATOM_MAP = AUDIT / "receptor_comparison" / "receptor_atom_comparison.csv"
RESIDUE_MAP = AUDIT / "receptor_comparison" / "receptor_residue_comparison.csv"
LIGAND_AUDIT = AUDIT / "ligand_protonation" / "ligand_protonation_audit.csv"
BOXES = {
    "7WC8": ROOT / "config" / "boxes" / "7wc8_box.json",
    "9LL8": ROOT / "config" / "boxes" / "9ll8_box.json",
    "9LLB": ROOT / "config" / "boxes" / "9llb_box.json",
}
RECEPTORS = {
    "7WC8": {
        "pdb": ROOT / "structures" / "7WC8_chainA_protein.pdb",
        "pdbqt": ROOT / "structures" / "7WC8_chainA_protein.pdbqt",
        "chain": "A",
        "docking_dir": ROOT / "results" / "7WC8_chainA_protein_functional",
    },
    "9LL8": {
        "pdb": ROOT / "structures" / "9LL8_chainR_protein.pdb",
        "pdbqt": ROOT / "structures" / "9LL8_chainR_protein.pdbqt",
        "chain": "R",
        "docking_dir": ROOT / "results" / "9LL8_chainR_protein_functional",
    },
    "9LLB": {
        "pdb": ROOT / "structures" / "9LLB_chainR_protein.pdb",
        "pdbqt": ROOT / "structures" / "9LLB_chainR_protein.pdbqt",
        "chain": "R",
        "docking_dir": ROOT / "results" / "9LLB_chainR_protein_functional",
    },
}
PDBQT_INPUTS = {
    "7WC8": ROOT / "structures" / "7WC8_chainA_protein.pdbqt",
    "9LL8": ROOT / "structures" / "9LL8_chainR_protein.pdbqt",
    "9LLB": ROOT / "structures" / "9LLB_chainR_protein.pdbqt",
}
BOX_PARAMETERS = {"exhaustiveness": 8, "num_modes": 9, "energy_range": 3, "cpu": 4}
PILOT_LIGANDS_PER_CLASS = 6
PH = 7.4
AMBER_ATOM_TYPES = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "S": "S", "SA": "S", "P": "P",
    "F": "F", "Cl": "Cl", "CL": "Cl", "Br": "Br", "BR": "Br",
    "I": "I", "HD": "H", "HS": "H",
}
BACKBONE = {"N", "CA", "C", "O", "CB"}
WT_REPAIRS = {
    ("7WC8", "162"): "SER",
    ("7WC8", "164"): "MET",
    ("7WC8", "372"): "SER",
}
# Joint template screening showed that this graft pair clashes when both
# sidechains come from 7WC8; using 9LLB for Lys146 leaves Asp218 from 7WC8
# and removes the close NZ-OD2 contact while retaining local backbone fits.
PREFERRED_TEMPLATE_DONORS = {("9LL8", "146"): "9LLB"}

sys.path.insert(0, str(ROOT / "scripts"))
import audit_docking_inputs as audit  # noqa: E402
import run_prolif_corrected as prolif_pipeline  # noqa: E402


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
    *,
    lineterminator: str = "\r\n",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Refusing to write an empty CSV: {path}")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
            extrasaction="ignore",
            lineterminator=lineterminator,
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes"}:
            return True
        if normalized in {"false", "0", "no", ""}:
            return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def raw_residue_key(line: str) -> tuple[str, str, str, str]:
    return (
        line[21:22].strip(),
        line[22:26].strip(),
        line[26:27].strip(),
        line[17:20].strip(),
    )


def atom_type(line: str) -> str:
    types = line[77:].split()
    if not types:
        raise ValueError(f"Missing AutoDock atom type: {line}")
    return types[0]


def atom_record(line: str) -> dict[str, Any]:
    atype = atom_type(line)
    element = AMBER_ATOM_TYPES.get(atype)
    if element is None:
        raise ValueError(f"Unsupported PDBQT atom type {atype!r}")
    return {
        "line": line,
        "atom_name": line[12:16].strip(),
        "element": element,
        "atom_type": atype,
        "residue_key": raw_residue_key(line),
        "xyz": np.asarray(
            [float(line[start:start + 8]) for start in (30, 38, 46)],
            dtype=float,
        ),
        "charge": float(line[70:76]),
    }


def format_pdbqt_atom(
    atom: dict[str, Any],
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resid: int,
    icode: str,
) -> str:
    x, y, z = (float(value) for value in atom["xyz"])
    return (
        f"ATOM  {serial:5d} {name:>4.4s} {resname:>3.3s} {chain[:1]:1s}"
        f"{resid:4d}{icode[:1]:1s}   {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.0:6.2f}{0.0:6.2f}    {float(atom['charge']):6.3f} {atom['atom_type']:<2s}"
    )


def format_pdb_atom(
    atom: dict[str, Any],
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resid: int,
    icode: str,
) -> str:
    x, y, z = (float(value) for value in atom["xyz"])
    atom_name = atom["name"][:4]
    atom_field = (
        f" {atom_name:<3s}"
        if len(atom["element"].strip()) == 1 and not atom_name[:1].isdigit()
        else f"{atom_name:<4s}"
    )
    return (
        f"ATOM  {serial:5d} {atom_field} {resname:>3.3s} {chain[:1]:1s}"
        f"{resid:4d}{icode[:1]:1s}   {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.0:6.2f}{0.0:6.2f}          {atom['element']:>2s}  "
    )


def kabsch_transform(
    source: np.ndarray, target: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] != 3:
        raise ValueError("Kabsch coordinate arrays must have matching Nx3 shapes")
    if len(source) < 3:
        raise ValueError("At least three matched atoms are required for alignment")
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    u, _, vt = np.linalg.svd((source - source_center).T @ (target - target_center))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    translation = target_center - source_center @ rotation
    fitted = source @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def residue_and_atom_maps(
    receptor: str,
    residue_rows: list[dict[str, str]],
    atom_rows: list[dict[str, str]],
) -> tuple[dict[tuple[str, str, str, str], dict[str, str]], dict]:
    residue_map = {}
    for row in residue_rows:
        if row["receptor"] != receptor:
            continue
        key = (
            row["pdbqt_chain"],
            row["pdbqt_residue"],
            row.get("pdbqt_insertion_code", ""),
            row["pdbqt_resname"],
        )
        if key in residue_map:
            raise ValueError(f"Ambiguous PDBQT residue mapping for {receptor}: {key}")
        residue_map[key] = row
    atom_map = {}
    for row in atom_rows:
        if row["receptor"] != receptor:
            continue
        key = (
            row["pdbqt_chain"],
            row["pdbqt_residue"],
            row.get("pdbqt_insertion_code", ""),
            row["pdbqt_resname"],
            row["pdbqt_atom"],
        )
        if key in atom_map:
            raise ValueError(f"Ambiguous PDBQT atom mapping for {receptor}: {key}")
        atom_map[key] = row["source_atom"]
    return residue_map, atom_map


def map_pdbqt_atoms(
    receptor: str,
    path: Path,
    residue_map: dict[tuple[str, str, str, str], dict[str, str]],
    atom_map: dict,
) -> tuple[list[dict[str, Any]], dict[tuple[str, str, str, str], list[dict[str, Any]]]]:
    records = []
    groups: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = atom_record(line)
        raw_key = atom["residue_key"]
        row = residue_map.get(raw_key)
        if row is None:
            raise ValueError(f"{receptor} PDBQT residue lacks source mapping: {raw_key}")
        atom["source_row"] = row
        atom["source_key"] = (
            row["source_chain"],
            row["source_residue"],
            row["source_insertion_code"],
            row["source_resname"],
        )
        atom["source_name"] = atom_map.get((*raw_key, atom["atom_name"]), atom["atom_name"])
        records.append(atom)
        groups[raw_key].append(atom)
    return records, groups


def canonical_groups(
    residue_rows: list[dict[str, str]],
    raw_groups: dict[tuple[str, str, str, str], list[dict[str, Any]]],
    receptor: str,
) -> dict[int, tuple[dict[str, str], dict[str, dict[str, Any]]]]:
    result = {}
    for row in residue_rows:
        if row["receptor"] != receptor or not row["source_uniprot_position"]:
            continue
        raw_key = (
            row["pdbqt_chain"],
            row["pdbqt_residue"],
            row.get("pdbqt_insertion_code", ""),
            row["pdbqt_resname"],
        )
        position = int(row["source_uniprot_position"])
        atoms = {}
        for atom in raw_groups.get(raw_key, []):
            if atom["element"] != "H":
                atoms[atom["source_name"]] = atom
        if position in result:
            continue
        result[position] = (row, atoms)
    return result


def local_alignment(
    target_receptor: str,
    target_position: int,
    donor_receptor: str,
    target_canonical: dict,
    donor_canonical: dict,
) -> tuple[np.ndarray, np.ndarray, float, int]:
    source_points = []
    target_points = []
    if target_position in target_canonical and target_position in donor_canonical:
        _, target_atoms = target_canonical[target_position]
        _, donor_atoms = donor_canonical[target_position]
        for name in ("N", "CA", "C", "O", "CB"):
            if name in donor_atoms and name in target_atoms:
                source_points.append(donor_atoms[name]["xyz"])
                target_points.append(target_atoms[name]["xyz"])
    if len(source_points) < 3:
        source_points = []
        target_points = []
        for position in range(target_position - 6, target_position + 7):
            if position not in target_canonical or position not in donor_canonical:
                continue
            _, target_atoms = target_canonical[position]
            _, donor_atoms = donor_canonical[position]
            for name in ("N", "CA", "C", "O", "CB"):
                if name in donor_atoms and name in target_atoms:
                    source_points.append(donor_atoms[name]["xyz"])
                    target_points.append(target_atoms[name]["xyz"])
    if len(source_points) < 3:
        raise ValueError(
            f"Too few residue-local backbone/CB anchors for {donor_receptor}->{target_receptor} "
            f"at {target_position}: {len(source_points)}"
        )
    rotation, translation, rmsd = kabsch_transform(
        np.asarray(source_points), np.asarray(target_points)
    )
    return rotation, translation, rmsd, len(source_points)


def sidechain_hydrogens(
    donor_atoms: list[dict[str, Any]], sidechain_names: set[str]
) -> list[dict[str, Any]]:
    heavy = [atom for atom in donor_atoms if atom["element"] != "H"]
    result = []
    for atom in donor_atoms:
        if atom["element"] != "H":
            continue
        nearest = min(
            heavy,
            key=lambda item: float(np.linalg.norm(atom["xyz"] - item["xyz"])),
            default=None,
        )
        if (
            nearest is not None
            and nearest["source_name"] in sidechain_names
            and np.linalg.norm(atom["xyz"] - nearest["xyz"]) <= 1.35
        ):
            result.append(atom)
    return result


def make_receptor(
    receptor: str,
    residue_rows: list[dict[str, str]],
    atom_rows: list[dict[str, str]],
    all_residue_rows: list[dict[str, str]],
    all_atom_rows: list[dict[str, str]],
    output_dir: Path,
) -> tuple[dict[str, Any], Path, Path]:
    target_residue_map, target_atom_map = residue_and_atom_maps(
        receptor, residue_rows, atom_rows
    )
    raw_records, raw_groups = map_pdbqt_atoms(
        receptor, PDBQT_INPUTS[receptor], target_residue_map, target_atom_map
    )
    target_canonical = canonical_groups(residue_rows, raw_groups, receptor)

    donors: dict[str, tuple[dict, list[dict], dict]] = {}
    for donor_name in RECEPTORS:
        donor_residue_map, donor_atom_map = residue_and_atom_maps(
            donor_name, all_residue_rows, all_atom_rows
        )
        donor_records, donor_groups = map_pdbqt_atoms(
            donor_name, PDBQT_INPUTS[donor_name], donor_residue_map, donor_atom_map
        )
        donors[donor_name] = (
            donor_residue_map,
            donor_records,
            donor_groups,
        )

    normalized = []
    for atom in raw_records:
        row = atom["source_row"]
        normalized.append(
            {
                **atom,
                "name": atom["source_name"] if atom["element"] != "H" else atom["atom_name"],
                "resname": row["source_resname"],
                "chain": row["source_chain"],
                "resid": int(row["source_residue"]),
                "icode": row["source_insertion_code"],
                "canonical": row["source_uniprot_position"],
                "original_source_name": row["source_resname"],
            }
        )

    additions: list[dict[str, Any]] = []
    removals: set[int] = set()
    repair_rows = []
    target_rows = [
        row for row in residue_rows
        if row["receptor"] == receptor
        and row["source_uniprot_position"]
        and float(row["minimum_heavy_atom_distance_to_box_angstrom"]) <= 6.0
    ]
    mutation_by_position = {
        position: resname
        for (rec, position), resname in WT_REPAIRS.items()
        if rec == receptor
    }
    work = []
    for row in target_rows:
        position = int(row["source_uniprot_position"])
        target_name = mutation_by_position.get(str(position))
        missing = {
            atom for atom in row["expected_heavy_atoms_missing_in_source"].split(";")
            if atom
        }
        if target_name or missing:
            work.append((row, target_name, missing))
    for row, desired_resname, missing in work:
        position = int(row["source_uniprot_position"])
        original_resname = row["source_resname"]
        desired_resname = desired_resname or original_resname
        donor_candidates = []
        for donor_name in RECEPTORS:
            if donor_name == receptor:
                continue
            donor_match = next(
                (
                    candidate for candidate in all_residue_rows
                    if candidate["receptor"] == donor_name
                    and candidate["source_uniprot_position"] == str(position)
                    and candidate["source_resname"] == desired_resname
                    and not candidate["expected_heavy_atoms_missing_in_source"]
                ),
                None,
            )
            if donor_match is None:
                continue
            donor_key = (
                donor_match["pdbqt_chain"],
                donor_match["pdbqt_residue"],
                donor_match.get("pdbqt_insertion_code", ""),
                donor_match["pdbqt_resname"],
            )
            donor_atoms = donors[donor_name][2].get(donor_key, [])
            donor_heavy_names = {
                atom["source_name"] for atom in donor_atoms if atom["element"] != "H"
            }
            required_sidechain = set(audit.EXPECTED_HEAVY[desired_resname]) - BACKBONE | {"CB"}
            if not required_sidechain.issubset(donor_heavy_names):
                continue
            try:
                rotation, translation, rmsd, anchor_count = local_alignment(
                    receptor, position, donor_name, target_canonical,
                    canonical_groups(
                        all_residue_rows, donors[donor_name][2], donor_name
                    ),
                )
            except ValueError:
                continue
            donor_candidates.append(
                (
                    rmsd, donor_name, donor_match, donor_atoms,
                    rotation, translation, anchor_count,
                )
            )
        if not donor_candidates:
            repair_rows.append(
                {
                    "receptor": receptor,
                    "uniprot_position": position,
                    "original_resname": original_resname,
                    "prepared_resname": desired_resname,
                    "issue": "no_complete_experimental_template",
                    "atoms_added": "",
                    "donor_receptor": "",
                    "local_backbone_cb_rmsd_angstrom": "",
                    "minimum_new_atom_clearance_angstrom": "",
                    "status": "not_reconstructed",
                }
            )
            continue

        target_key = (
            row["source_chain"], row["source_residue"],
            row["source_insertion_code"], original_resname,
        )
        target_residue_atoms = [
            atom for atom in normalized
            if (
                atom["chain"], str(atom["resid"]), atom["icode"], atom["resname"]
            ) == target_key
        ]
        existing_names = {
            atom["name"] for atom in target_residue_atoms
            if atom["element"] != "H"
        }
        replace_mutation = desired_resname != original_resname
        if replace_mutation:
            remove_names = {
                atom["name"] for atom in target_residue_atoms
                if atom["element"] != "H" and atom["name"] not in BACKBONE
            }
            keep_heavy = [
                atom for atom in target_residue_atoms
                if atom["element"] != "H" and atom["name"] in BACKBONE
            ]
            for atom in target_residue_atoms:
                if atom["element"] == "H":
                    nearest = min(
                        (
                            (float(np.linalg.norm(atom["xyz"] - heavy["xyz"])), heavy)
                            for heavy in target_residue_atoms
                            if heavy["element"] != "H"
                        ),
                        key=lambda item: item[0],
                        default=(math.inf, None),
                    )[1]
                    if nearest and nearest["name"] not in BACKBONE:
                        remove_names.add(atom["name"])
            removed_line_indexes = {
                id(atom) for atom in target_residue_atoms
                if atom["name"] in remove_names
            }
            removals.update(
                id(atom) for atom in target_residue_atoms
                if id(atom) in removed_line_indexes
            )
            existing_names = {atom["name"] for atom in keep_heavy}
        required_sidechain = set(audit.EXPECTED_HEAVY[desired_resname]) - BACKBONE
        needed_names = required_sidechain - existing_names
        chosen = None
        preferred_donor = PREFERRED_TEMPLATE_DONORS.get(
            (receptor, str(position))
        )
        for candidate in sorted(
            donor_candidates,
            key=lambda item: (
                item[1] != preferred_donor if preferred_donor else False,
                item[0],
                item[1],
            ),
        ):
            rmsd, donor_name, donor_match, donor_atoms, rotation, translation, anchors = candidate
            donor_heavy = [
                atom for atom in donor_atoms if atom["element"] != "H"
                and atom["source_name"] in needed_names
            ]
            transformed = []
            for atom in donor_heavy:
                transformed.append(
                    {
                        **atom,
                        "xyz": atom["xyz"] @ rotation + translation,
                        "name": atom["source_name"],
                        "resname": desired_resname,
                        "chain": row["source_chain"],
                        "resid": int(row["source_residue"]),
                        "icode": row["source_insertion_code"],
                        "canonical": str(position),
                    }
                )
            external_heavy = [
                atom for atom in normalized
                if atom["element"] != "H"
                and atom["canonical"] != str(position)
                and id(atom) not in removals
            ] + additions
            clearances = [
                float(np.linalg.norm(atom["xyz"] - other["xyz"]))
                for atom in transformed for other in external_heavy
            ]
            clearance = min(clearances, default=math.inf)
            if clearance >= 2.0:
                chosen = (*candidate, transformed, clearance)
                break
            if chosen is None:
                chosen = (*candidate, transformed, clearance)
        if chosen is None:
            raise RuntimeError(f"No structural graft candidate for {receptor}:{position}")
        rmsd, donor_name, donor_match, donor_atoms, rotation, translation, anchors, transformed, clearance = chosen
        if len(transformed) != len(needed_names):
            raise ValueError(
                f"Incomplete graft at {receptor}:{position}; "
                f"needed={sorted(needed_names)}, obtained={sorted(atom['name'] for atom in transformed)}"
            )
        new_atom_names = {atom["name"] for atom in transformed}
        for atom in sidechain_hydrogens(donor_atoms, new_atom_names):
            additions.append(
                {
                    **atom,
                    "xyz": atom["xyz"] @ rotation + translation,
                    "name": atom["atom_name"],
                    "resname": desired_resname,
                    "chain": row["source_chain"],
                    "resid": int(row["source_residue"]),
                    "icode": row["source_insertion_code"],
                    "canonical": str(position),
                }
            )
        additions.extend(transformed)
        for atom in normalized:
            if (
                atom["element"] != "H"
                and atom["canonical"] == str(position)
                and atom["name"] in BACKBONE
            ):
                atom["resname"] = desired_resname
        for atom in normalized:
            if atom["element"] == "H" and atom["canonical"] == str(position):
                atom["resname"] = desired_resname
        repair_rows.append(
            {
                "receptor": receptor,
                "uniprot_position": position,
                "original_resname": original_resname,
                "prepared_resname": desired_resname,
                "issue": "engineered_mutation" if replace_mutation else "incomplete_sidechain",
                "template_selection_reason": (
                    "joint sidechain fit with ASP218 minimizes a close contact"
                    if donor_name == preferred_donor else
                    "lowest local fit with acceptable steric clearance"
                ),
                "atoms_added": ";".join(sorted(new_atom_names)),
                "donor_receptor": donor_name,
                "local_backbone_cb_rmsd_angstrom": round(rmsd, 3),
                "local_alignment_atoms": anchors,
                "minimum_new_atom_clearance_angstrom": (
                    round(clearance, 3) if math.isfinite(clearance) else ""
                ),
                "status": (
                    "grafted_low_confidence" if rmsd > 1.0 or clearance < 2.0
                    else "grafted_from_local_experimental_template"
                ),
            }
        )

    base = [atom for atom in normalized if id(atom) not in removals]
    atom_groups: dict[tuple[str, int, str, str], list[dict[str, Any]]] = {}
    residue_order = []
    for atom in base + additions:
        key = (atom["chain"], atom["resid"], atom["icode"], atom["resname"])
        if key not in atom_groups:
            atom_groups[key] = []
            residue_order.append(key)
        atom_groups[key].append(atom)
    all_atoms = [
        atom
        for key in residue_order
        for atom in atom_groups[key]
    ]
    residues = defaultdict(set)
    for atom in all_atoms:
        if atom["element"] != "H":
            residues[(atom["chain"], atom["resid"], atom["icode"], atom["resname"])].add(atom["name"])
    for repair in repair_rows:
        if repair["status"].startswith("grafted"):
            position = str(repair["uniprot_position"])
            row = next(
                item for item in residue_rows
                if item["receptor"] == receptor
                and item["source_uniprot_position"] == position
            )
            key = (
                row["source_chain"], int(row["source_residue"]),
                row["source_insertion_code"], repair["prepared_resname"],
            )
            missing_after = set(audit.EXPECTED_HEAVY[repair["prepared_resname"]]) - residues[key]
            if missing_after:
                raise ValueError(
                    f"Corrected {receptor}:{position} remains incomplete: {sorted(missing_after)}"
                )
            grafted_names = set(repair["atoms_added"].split(";"))
            grafted_atoms = [
                atom for atom in all_atoms
                if atom["canonical"] == position
                and atom["name"] in grafted_names
                and atom["element"] != "H"
            ]
            external_atoms = [
                atom for atom in all_atoms
                if atom["canonical"] != position
                and atom["element"] != "H"
            ]
            closest = min(
                (
                    (
                        float(np.linalg.norm(atom["xyz"] - other["xyz"])),
                        atom,
                        other,
                    )
                    for atom in grafted_atoms
                    for other in external_atoms
                ),
                key=lambda item: item[0],
                default=None,
            )
            if closest is not None:
                clearance, grafted_atom, external_atom = closest
                repair["minimum_new_atom_clearance_angstrom"] = round(clearance, 3)
                repair["closest_external_atom"] = (
                    f"{grafted_atom['name']}:{external_atom['resname']}"
                    f"{external_atom['canonical']}:{external_atom['name']}"
                )
                if clearance < 2.0:
                    repair["status"] = "grafted_low_confidence"
    seen = set()
    for atom in all_atoms:
        if atom["element"] == "H":
            continue
        key = (atom["chain"], atom["resid"], atom["icode"], atom["resname"], atom["name"])
        if key in seen:
            raise ValueError(f"Duplicate heavy atom in corrected receptor {receptor}: {key}")
        seen.add(key)

    pdbqt_path = output_dir / "receptor_preparation" / f"{receptor}_recovered.pdbqt"
    pdb_path = output_dir / "receptor_preparation" / f"{receptor}_recovered.pdb"
    pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        f"REMARK  Localized pH 7.4 receptor recovery: {receptor}",
        "REMARK  Standard Vina scoring; atom types/charges from experimental receptor PDBQT templates",
    ]
    pdbqt_lines = header + [
        format_pdbqt_atom(
            atom, index, atom["name"], atom["resname"], atom["chain"],
            int(atom["resid"]), atom["icode"],
        )
        for index, atom in enumerate(all_atoms, start=1)
    ] + ["END"]
    pdbqt_path.write_text("\n".join(pdbqt_lines) + "\n", encoding="utf-8")
    heavy = [atom for atom in all_atoms if atom["element"] != "H"]
    pdb_lines = [
        format_pdb_atom(
            atom, index, atom["name"], atom["resname"], atom["chain"],
            int(atom["resid"]), atom["icode"],
        )
        for index, atom in enumerate(heavy, start=1)
    ] + ["END"]
    pdb_path.write_text("\n".join(pdb_lines) + "\n", encoding="utf-8")

    unresolved = [
        row for row in repair_rows if row["status"] == "not_reconstructed"
    ]
    unresolved_ids = {
        (row["receptor"], str(row["uniprot_position"])) for row in unresolved
    }
    unresolved_rows = [
        row for row in target_rows
        if row["expected_heavy_atoms_missing_in_source"]
        and (receptor, row["source_uniprot_position"]) in unresolved_ids
    ]
    expected_unresolved = [
        row for row in unresolved_rows
        if (
            float(row["minimum_distance_to_pilot_ligand_angstrom"] or math.inf) <= 8.0
            or float(row["minimum_heavy_atom_distance_to_box_angstrom"]) == 0.0
        )
    ]
    deferred = [
        {
            "uniprot_position": row["source_uniprot_position"],
            "resname": row["source_resname"],
            "missing_atoms": row["expected_heavy_atoms_missing_in_source"],
            "distance_to_box_angstrom": float(
                row["minimum_heavy_atom_distance_to_box_angstrom"]
            ),
            "distance_to_nearest_pilot_pose_angstrom": float(
                row["minimum_distance_to_pilot_ligand_angstrom"] or math.inf
            ),
        }
        for row in unresolved_rows if row not in expected_unresolved
    ]
    return {
        "receptor": receptor,
        "source_pdb": str(RECEPTORS[receptor]["pdb"].relative_to(ROOT)),
        "source_pdbqt": str(PDBQT_INPUTS[receptor].relative_to(ROOT)),
        "prepared_pdb": str(pdb_path.relative_to(ROOT)),
        "prepared_pdbqt": str(pdbqt_path.relative_to(ROOT)),
        "source_sha256": {
            "pdb": sha256(RECEPTORS[receptor]["pdb"]),
            "pdbqt": sha256(PDBQT_INPUTS[receptor]),
        },
        "prepared_sha256": {
            "pdb": sha256(pdb_path),
            "pdbqt": sha256(pdbqt_path),
        },
        "heavy_atoms_before": sum(atom["element"] != "H" for atom in normalized),
        "heavy_atoms_after": len(heavy),
        "heavy_atom_coordinate_changes": (
            sum(row["issue"] == "engineered_mutation" for row in repair_rows)
        ),
        "repairs": repair_rows,
        "unresolved": expected_unresolved,
        "deferred_peripheral_residues": deferred,
        "pdbqt_has_nonzero_partial_charges": any(
            abs(atom["charge"]) > 1e-6 for atom in all_atoms
        ),
    }, pdb_path, pdbqt_path


def add_protomer_coordinates(
    reference: Chem.Mol,
    protomer: Chem.Mol,
    mapping: dict[int, int],
) -> Chem.Mol:
    reference_heavy = Chem.RemoveHs(Chem.Mol(reference))
    source_conf = reference_heavy.GetConformer()
    conformer = Chem.Conformer(protomer.GetNumAtoms())
    for source_index, state_index in mapping.items():
        point = source_conf.GetAtomPosition(source_index)
        conformer.SetAtomPosition(state_index, point)
    for source_atom in reference_heavy.GetAtoms():
        if source_atom.GetSymbol() != "H" or source_atom.GetIsotope() == 0:
            continue
        source_parent = source_atom.GetNeighbors()[0].GetIdx()
        state_parent = mapping[source_parent]
        state_deuterium = next(
            (
                atom.GetIdx() for atom in protomer.GetAtomWithIdx(state_parent).GetNeighbors()
                if atom.GetSymbol() == "H" and atom.GetIsotope() == source_atom.GetIsotope()
            ),
            None,
        )
        if state_deuterium is None:
            raise ValueError("Isotopic hydrogen was not retained during protomer selection")
        conformer.SetAtomPosition(
            state_deuterium, source_conf.GetAtomPosition(source_atom.GetIdx())
        )
    result = Chem.Mol(protomer)
    result.RemoveAllConformers()
    result.AddConformer(conformer, assignId=True)
    result = Chem.AddHs(result, addCoords=True)
    Chem.SanitizeMol(result)
    return result


def map_pdbqt_to_protomer(
    reference: Chem.Mol,
    protomer: Chem.Mol,
    reference_to_protomer: dict[int, int],
    prepared_pdbqt: Path,
) -> tuple[dict[int, int], list[int], dict[str, Any]]:
    heavy_reference = {
        atom.GetIdx() for atom in reference.GetAtoms()
        if atom.GetSymbol() != "H"
    }
    heavy_protomer = {
        atom.GetIdx() for atom in protomer.GetAtoms()
        if atom.GetSymbol() != "H"
    }
    mapped_heavy = {
        reference_index: protomer_index
        for reference_index, protomer_index in reference_to_protomer.items()
        if reference.GetAtomWithIdx(reference_index).GetSymbol() != "H"
    }
    if set(mapped_heavy) != heavy_reference or set(mapped_heavy.values()) != heavy_protomer:
        raise ValueError("SDF-to-protomer mapping does not bijectively cover heavy atoms")
    for reference_index, protomer_index in mapped_heavy.items():
        if (
            reference.GetAtomWithIdx(reference_index).GetAtomicNum()
            != protomer.GetAtomWithIdx(protomer_index).GetAtomicNum()
        ):
            raise ValueError("SDF-to-protomer mapping changes a heavy-atom element")

    protomer_coordinates = protomer.GetConformer().GetPositions()
    available = set(heavy_protomer)
    pdbqt_heavy = [
        atom for atom in prolif_pipeline.pdbqt_atoms(prepared_pdbqt)
        if atom["element"] != "H"
    ]
    ordered_indices = []
    maximum_delta = 0.0
    for record in pdbqt_heavy:
        coordinate = np.asarray([record[axis] for axis in ("x", "y", "z")])
        candidates = sorted(
            (
                float(np.linalg.norm(protomer_coordinates[index] - coordinate)),
                index,
            )
            for index in available
            if protomer.GetAtomWithIdx(index).GetSymbol() == record["element"]
        )
        if not candidates or candidates[0][0] > 0.002:
            raise ValueError(
                f"No coordinate-identical SDF atom for PDBQT serial {record['serial']}"
            )
        if len(candidates) > 1 and candidates[1][0] - candidates[0][0] < 0.0005:
            raise ValueError(
                f"Ambiguous coordinate mapping for PDBQT serial {record['serial']}"
            )
        delta, atom_index = candidates[0]
        available.remove(atom_index)
        ordered_indices.append(atom_index)
        maximum_delta = max(maximum_delta, delta)
    if available or len(ordered_indices) != len(heavy_protomer):
        raise ValueError("PDBQT heavy-atom set differs from selected protomer")
    return mapped_heavy, ordered_indices, {
        "reference_heavy_atoms": len(heavy_reference),
        "pdbqt_heavy_atoms": len(pdbqt_heavy),
        "maximum_coordinate_delta_angstrom": round(maximum_delta, 6),
        "mapping_status": "unique_by_element_and_heavy_atom_coordinates",
    }


def prepare_ligands(
    output_dir: Path,
    obabel: str,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    from run_prolif_corrected import read_sdf_records, select_protomer

    records = read_sdf_records(LIGAND_SDF)
    selection_rows = {row["ligand_id"]: row for row in read_csv(SELECTION)}
    audit_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_csv(LIGAND_AUDIT):
        audit_rows[row["ligand_id"]].append(row)
    out_ligands = output_dir / "ligand_preparation" / "pdbqt"
    tmp_sdf = output_dir / "ligand_preparation" / "sdf"
    out_ligands.mkdir(parents=True, exist_ok=True)
    tmp_sdf.mkdir(parents=True, exist_ok=True)
    summary = []
    prepared = {}
    for ligand_id, reference in sorted(records.items()):
        selection = selection_rows.get(ligand_id, {})
        qc = next(
            (row for row in audit_rows.get(ligand_id, []) if row["microstate_index"] == "1"),
            {},
        )
        output_pdbqt = out_ligands / f"{ligand_id}.pdbqt"
        source_pdbqt = ROOT / "structures" / "ligands" / "pdbqt" / f"{ligand_id}.pdbqt"
        status = "prepared"
        error = ""
        report: dict[str, Any] = {}
        try:
            if not source_pdbqt.exists():
                raise FileNotFoundError("no original docking PDBQT; excluded by input QC")
            protomer, mapping, report = select_protomer(reference, PH)
            mol_3d = add_protomer_coordinates(reference, protomer, mapping)
            for key, value in reference.GetPropsAsDict().items():
                mol_3d.SetProp(str(key), str(value))
            sdf_path = tmp_sdf / f"{ligand_id}.sdf"
            Chem.MolToMolFile(mol_3d, str(sdf_path), forceV3000=True)
            result = subprocess.run(
                [obabel, "-isdf", str(sdf_path), "-opdbqt", "-O", str(output_pdbqt)],
                capture_output=True,
                text=True,
                timeout=60,
            )
            (output_dir / "logs" / f"ligand_{ligand_id}.log").write_text(
                result.stdout + result.stderr, encoding="utf-8"
            )
            if result.returncode != 0 or not output_pdbqt.exists():
                raise RuntimeError(f"Open Babel conversion failed: {result.stderr[-500:]}")
            atom_map, prepared_indices, map_report = map_pdbqt_to_protomer(
                reference, mol_3d, mapping, output_pdbqt
            )
            prepared[ligand_id] = {
                "reference": reference,
                "protomer": protomer,
                "protomer_mapping": mapping,
                "pdbqt": output_pdbqt,
                "converted": protomer,
                "atom_mapping": atom_map,
                "map_report": map_report,
                "prepared_indices": prepared_indices,
                "protonation_report": report,
            }
            if report["state_ambiguity"]:
                status = "prepared_state_requires_manual_review"
        except Exception as exc:
            status = "unavailable_or_manual_review"
            error = f"{type(exc).__name__}: {exc}"
            output_pdbqt.unlink(missing_ok=True)
        original_charge = int(selection.get("structure_charge", "0") or "0")
        new_charge = report.get("formal_charge")
        original_smiles = selection.get("structure_smiles", "")
        protonation_changed = (
            new_charge is not None and int(new_charge) != original_charge
        )
        summary.append(
            {
                "ligand_id": ligand_id,
                "pharmacological_label": selection.get("pharmacological_label", ""),
                "original_state_smiles": original_smiles,
                "selected_state_smiles": report.get("chosen_smiles", ""),
                "original_formal_charge": original_charge,
                "selected_formal_charge": new_charge if new_charge is not None else "",
                "basic_site_count": report.get("basic_site_count", ""),
                "carboxyl_oxygen_count": report.get("carboxyl_oxygen_count", ""),
                "candidate_microstate_count": report.get("candidate_states", ""),
                "state_ambiguity": report.get("state_ambiguity", ""),
                "candidate_states_json": json.dumps(
                    report.get("candidate_smiles", []), ensure_ascii=False
                ),
                "protonation_selection_rule": report.get("selection_policy", ""),
                "population_estimated": False,
                "pdbqt_available": output_pdbqt.exists(),
                "heavy_atom_mapping_status": (
                    prepared.get(ligand_id, {}).get("map_report", {}).get("mapping_status", "")
                ),
                "original_docking_state_differs_by_charge": protonation_changed,
                "ligand_decision": (
                    "REVISIÓN MANUAL"
                    if status != "prepared" else
                    "REPETIR"
                    if protonation_changed else "CONSERVAR"
                ),
                "preparation_status": status,
                "error": error,
            }
        )
    write_csv(output_dir / "ligand_preparation_summary.csv", summary)
    return summary, prepared


def select_pilot(
    prepared: dict[str, dict[str, Any]],
    selection_rows: list[dict[str, str]],
    ligand_summary: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    selection = {row["ligand_id"]: row for row in selection_rows}
    prep_summary = {row["ligand_id"]: row for row in ligand_summary}
    candidates: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for ligand_id, data in prepared.items():
        row = selection.get(ligand_id, {})
        label = row.get("pharmacological_label", "")
        if label not in {"agonist", "antagonist"}:
            continue
        molecule = Chem.RemoveHs(Chem.Mol(data["reference"]))
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)
        if not scaffold:
            scaffold = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        ring_count = sum(atom.GetIsAromatic() for atom in molecule.GetAtoms())
        score = (
            4 * bool(data["protonation_report"]["basic_site_count"])
            + 2 * bool(Chem.rdMolDescriptors.CalcTPSA(molecule) > 20)
            + bool(ring_count)
            + bool(data["protonation_report"]["candidate_states"] > 1)
            - 3 * bool(data["protonation_report"]["state_ambiguity"])
        )
        candidates[label].append(
            {
                "ligand_id": ligand_id,
                "label": label,
                "scaffold": scaffold,
                "score": score,
                "state": prep_summary[ligand_id],
                "basic_sites": data["protonation_report"]["basic_site_count"],
                "hbd": Chem.rdMolDescriptors.CalcNumHBD(molecule),
                "hba": Chem.rdMolDescriptors.CalcNumHBA(molecule),
                "aromatic_atoms": ring_count,
                "heavy_atoms": sum(atom.GetSymbol() != "H" for atom in molecule.GetAtoms()),
            }
        )
    pilot = []
    for label in ("agonist", "antagonist"):
        ordered = sorted(
            candidates[label],
            key=lambda row: (
                parse_bool(row["state"]["state_ambiguity"]),
                -row["score"],
                -row["heavy_atoms"],
                row["ligand_id"],
            ),
        )
        chosen = []
        used_scaffolds = set()
        for unique_pass in (True, False):
            for row in ordered:
                if row["ligand_id"] in {item["ligand_id"] for item in chosen}:
                    continue
                if unique_pass and row["scaffold"] in used_scaffolds:
                    continue
                chosen.append(row)
                used_scaffolds.add(row["scaffold"])
                if len(chosen) == PILOT_LIGANDS_PER_CLASS:
                    break
            if len(chosen) == PILOT_LIGANDS_PER_CLASS:
                break
        if len(chosen) < PILOT_LIGANDS_PER_CLASS:
            raise RuntimeError(
                f"Only {len(chosen)} prepared {label} ligands available for pilot"
            )
        pilot.extend(chosen)
    if len(pilot) < 12 or len(pilot) > 20:
        raise ValueError(f"Pilot size is outside the requested range: {len(pilot)}")
    return pilot


def run_pilot_docking(
    pilot: list[dict[str, Any]],
    prepared: dict[str, dict[str, Any]],
    output_dir: Path,
    vina: str,
) -> list[dict[str, Any]]:
    rows = []
    seed = 4512026
    for receptor in RECEPTORS:
        receptor_path = output_dir / "receptor_preparation" / f"{receptor}_recovered.pdbqt"
        box = json.loads(BOXES[receptor].read_text(encoding="utf-8"))
        for item in pilot:
            ligand_id = item["ligand_id"]
            ligand_path = prepared[ligand_id]["pdbqt"]
            run_dir = output_dir / "pilot_docking" / receptor
            for subdir in ("poses", "logs", "configs"):
                (run_dir / subdir).mkdir(parents=True, exist_ok=True)
            out_pose = run_dir / "poses" / f"{ligand_id}.pdbqt"
            log_path = run_dir / "logs" / f"{ligand_id}.log"
            config_path = run_dir / "configs" / f"{ligand_id}.txt"
            config = {
                "receptor": str(receptor_path.resolve()),
                "ligand": str(ligand_path.resolve()),
                "out": str(out_pose.resolve()),
                "center_x": box["center_x"],
                "center_y": box["center_y"],
                "center_z": box["center_z"],
                "size_x": box["size_x"],
                "size_y": box["size_y"],
                "size_z": box["size_z"],
                **BOX_PARAMETERS,
                "scoring": "vina",
                "seed": seed,
            }
            seed += 1
            config_path.write_text(
                "\n".join(f"{key} = {value}" for key, value in config.items()) + "\n",
                encoding="utf-8",
            )
            command = [
                vina,
                "--receptor", str(receptor_path),
                "--ligand", str(ligand_path),
                "--out", str(out_pose),
                "--center_x", str(box["center_x"]),
                "--center_y", str(box["center_y"]),
                "--center_z", str(box["center_z"]),
                "--size_x", str(box["size_x"]),
                "--size_y", str(box["size_y"]),
                "--size_z", str(box["size_z"]),
                "--exhaustiveness", str(BOX_PARAMETERS["exhaustiveness"]),
                "--num_modes", str(BOX_PARAMETERS["num_modes"]),
                "--energy_range", str(BOX_PARAMETERS["energy_range"]),
                "--cpu", str(BOX_PARAMETERS["cpu"]),
                "--seed", str(config["seed"]),
            ]
            result = subprocess.run(
                command, capture_output=True, text=True, timeout=1800
            )
            log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
            if result.returncode != 0 or not out_pose.exists():
                rows.append(
                    {
                        "receptor": receptor,
                        "ligand_id": ligand_id,
                        "pharmacological_label": item["label"],
                        "status": "failed",
                        "returncode": result.returncode,
                        "top_affinity_kcal_mol": "",
                        "pose_count": 0,
                        "seed": config["seed"],
                        "config": str(config_path.relative_to(ROOT)),
                        "pose": "",
                        "log": str(log_path.relative_to(ROOT)),
                        "error": (result.stderr or result.stdout)[-500:],
                    }
                )
                continue
            scores = re.findall(r"^\s*\d+\s+(-?\d+\.\d+)", result.stdout, re.MULTILINE)
            models = prolif_pipeline.parse_pose_models(out_pose)
            rows.append(
                {
                    "receptor": receptor,
                    "ligand_id": ligand_id,
                    "pharmacological_label": item["label"],
                    "status": "ok",
                    "returncode": result.returncode,
                    "top_affinity_kcal_mol": float(scores[0]) if scores else "",
                    "pose_count": len(models),
                    "seed": config["seed"],
                    "config": str(config_path.relative_to(ROOT)),
                    "pose": str(out_pose.relative_to(ROOT)),
                    "log": str(log_path.relative_to(ROOT)),
                    "error": "",
                }
            )
    write_csv(output_dir / "pilot_docking" / "pilot_docking_results.csv", rows)
    return rows


def process_pilot_prolif(
    pilot: list[dict[str, Any]],
    prepared: dict[str, dict[str, Any]],
    docking_rows: list[dict[str, Any]],
    output_dir: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    for receptor in RECEPTORS:
        prolif_pipeline.RECEPTORS[receptor]["pdb"] = (
            output_dir / "receptor_preparation" / f"{receptor}_recovered.pdb"
        )
        prolif_pipeline.RECEPTORS[receptor]["pdbqt"] = (
            output_dir / "receptor_preparation" / f"{receptor}_recovered.pdbqt"
        )
    proteins = {}
    receptor_reports = {}
    for receptor in RECEPTORS:
        proteins[receptor], receptor_reports[receptor] = prolif_pipeline.prepare_receptor(
            receptor, PH, output_dir / "pilot_prolif"
        )

    rows_by_receptor = defaultdict(list)
    docking_by_key = {
        (row["receptor"], row["ligand_id"]): row
        for row in docking_rows if row["status"] == "ok"
    }
    pose_records = []
    for receptor in RECEPTORS:
        for item in pilot:
            ligand_id = item["ligand_id"]
            docking = docking_by_key.get((receptor, ligand_id))
            if docking is None:
                continue
            state = prepared[ligand_id]
            converted = state["converted"]
            heavy_indices = state["prepared_indices"]
            atom_mapping = {
                "ref_to_prepared": state["atom_mapping"],
                "prepared_heavy_indices": heavy_indices,
                "elements": {
                    index: converted.GetAtomWithIdx(index).GetSymbol()
                    for index in heavy_indices
                },
            }
            pose_path = ROOT / docking["pose"]
            models = prolif_pipeline.parse_pose_models(pose_path)
            complex_record = {
                "ligand_id": ligand_id,
                "pharmacological_label": item["label"],
                "pose": str(pose_path),
            }
            for model in models:
                processed = prolif_pipeline.process_one_pose(
                    receptor,
                    complex_record,
                    model["model"],
                    model["atoms"],
                    state["reference"],
                    state["protomer"],
                    state["protomer_mapping"],
                    atom_mapping,
                    proteins[receptor],
                    output_dir / "pilot_prolif",
                )
                payload = processed["payload"]
                rows_by_receptor[receptor].append(
                    {
                        "pose_id": f"{receptor}_{ligand_id}_model_{model['model']:02d}",
                        "ligand_id": ligand_id,
                        "pharmacological_label": item["label"],
                        "pose_model": model["model"],
                        "interaction_event_counts": payload["interaction_event_counts"],
                        "active_features": sorted(processed["active_features"]),
                        "bits": processed["bits"],
                        "features": processed["features"],
                    }
                )
                pose_records.append(payload)

    reproducibility = {"tested": False, "identical": False}
    first = next(
        (
            (
                receptor,
                item,
                prepared[item["ligand_id"]],
                docking_by_key[(receptor, item["ligand_id"])],
            )
            for receptor in RECEPTORS for item in pilot
            if (receptor, item["ligand_id"]) in docking_by_key
        ),
        None,
    )
    if first:
        receptor, item, state, docking = first
        model = prolif_pipeline.parse_pose_models(ROOT / docking["pose"])[0]
        converted = state["converted"]
        indices = state["prepared_indices"]
        mapping = {
            "ref_to_prepared": state["atom_mapping"],
            "prepared_heavy_indices": indices,
            "elements": {
                index: converted.GetAtomWithIdx(index).GetSymbol()
                for index in indices
            },
        }
        complex_record = {
            "ligand_id": item["ligand_id"],
            "pharmacological_label": item["label"],
            "pose": docking["pose"],
        }
        repeated = prolif_pipeline.process_one_pose(
            receptor, complex_record, model["model"], model["atoms"],
            state["reference"], state["protomer"], state["protomer_mapping"],
            mapping, proteins[receptor], output_dir / "pilot_prolif" / "reproducibility",
        )
        first_row = rows_by_receptor[receptor][0]
        reproducibility = {
            "tested": True,
            "identical": bool(
                np.array_equal(first_row["bits"], repeated["bits"])
                and first_row["active_features"] == sorted(repeated["active_features"])
            ),
            "complex": f"{receptor}:{item['ligand_id']}:model{model['model']}",
        }

    summary = {
        "n_poses": len(pose_records),
        "reproducibility": reproducibility,
        "receptor_preparation": receptor_reports,
        "receptors": {},
        "interaction_categories": list(prolif_pipeline.INTERACTIONS),
    }
    for receptor, rows in rows_by_receptor.items():
        feature_names = sorted({feature for row in rows for feature in row["features"]})
        feature_index = {feature: index for index, feature in enumerate(feature_names)}
        matrix = np.zeros((len(rows), len(feature_names)), dtype=np.uint8)
        for row_index, row in enumerate(rows):
            for feature in row["active_features"]:
                matrix[row_index, feature_index[feature]] = 1
        counts = Counter()
        for row in rows:
            counts.update(row["interaction_event_counts"])
        matrix_dir = output_dir / "pilot_prolif" / "matrices"
        matrix_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            matrix_dir / f"IFP_{receptor}.npz",
            matrix=matrix,
            pose_ids=np.asarray([row["pose_id"] for row in rows], dtype=str),
            feature_names=np.asarray(feature_names, dtype=str),
            interaction_types=np.asarray(prolif_pipeline.INTERACTIONS, dtype=str),
        )
        csv_rows = []
        for index, row in enumerate(rows):
            csv_rows.append(
                {
                    "pose_id": row["pose_id"],
                    "ligand_id": row["ligand_id"],
                    "pharmacological_label": row["pharmacological_label"],
                    "pose_model": row["pose_model"],
                    **{
                        feature: int(matrix[index, feature_index[feature]])
                        for feature in feature_names
                    },
                }
            )
        write_csv(matrix_dir / f"IFP_{receptor}.csv", csv_rows)
        summary["receptors"][receptor] = {
            "pose_count": len(rows),
            "interaction_event_counts": {
                name: int(counts[name]) for name in prolif_pipeline.INTERACTIONS
            },
            "category_positive_pose_counts": {
                name: sum(
                    row["interaction_event_counts"].get(name, 0) > 0
                    for row in rows
                )
                for name in prolif_pipeline.INTERACTIONS
            },
            "diagnostic_ASP155_events": sum(
                1 for payload in pose_records
                if payload["receptor_key"] == receptor
                for event in payload["events"]
                if "ASP155" in event.get("protein_residue", "")
            ),
        }
    (output_dir / "pilot_prolif" / "pilot_prolif_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return pose_records, summary


def create_decisions(
    ligand_summary: list[dict[str, Any]],
    receptor_preparations: dict[str, dict[str, Any]],
    output_dir: Path,
) -> list[dict[str, Any]]:
    rows = []
    summary = {row["ligand_id"]: row for row in ligand_summary}
    for receptor, prepared in receptor_preparations.items():
        has_structural_change = any(
            row["status"].startswith("grafted")
            for row in prepared["repairs"]
        )
        source_dir = RECEPTORS[receptor]["docking_dir"] / "poses"
        for ligand_id, ligand in sorted(summary.items()):
            old_pose = source_dir / f"{ligand_id}.pdbqt"
            if not old_pose.exists():
                decision = "NO DISPONIBLE"
                problem = "No hay pose PDBQT original; ligando ausente o excluido"
                justification = "No se inventa una estructura ni un resultado histórico."
            elif ligand["preparation_status"] == "unavailable_or_manual_review":
                decision = "REVISIÓN MANUAL"
                problem = ligand["error"]
                justification = "La representación a pH 7.4 no se pudo validar."
            elif parse_bool(ligand["state_ambiguity"]):
                decision = "REVISIÓN MANUAL"
                problem = "Más de un estado de microestado plausible a pH 7.4"
                justification = (
                    "Se generó y registró un estado representativo para el piloto; "
                    "la elección requiere sensibilidad dirigida antes del redocking completo."
                )
            elif (
                has_structural_change
                or parse_bool(ligand["original_docking_state_differs_by_charge"])
            ):
                decision = "REPETIR"
                problems = []
                if has_structural_change:
                    problems.append("receptor del bolsillo preparado/corregido")
                if ligand["original_docking_state_differs_by_charge"]:
                    problems.append("carga del ligando cambia frente al estado original")
                problem = "; ".join(problems)
                justification = (
                    "El cambio del receptor es compartido por todos los ligandos de este "
                    "receptor; no depende de la pose final. Se mantienen los parámetros "
                    "Vina comparables y se asigna un seed explícito."
                )
            else:
                decision = "CONSERVAR"
                problem = "No se detectó cambio químico en receptor o carga del ligando"
                justification = "La comparación del estado y del input no exige repetir."
            rows.append(
                {
                    "ligand_id": ligand_id,
                    "receptor": receptor,
                    "problem": problem,
                    "decision": decision,
                    "justification": justification,
                    "original_pose_available": old_pose.exists(),
                    "original_charge": ligand["original_formal_charge"],
                    "selected_charge": ligand["selected_formal_charge"],
                    "selected_state_smiles": ligand["selected_state_smiles"],
                }
            )
    write_csv(output_dir / "redocking_decisions.csv", rows)
    return rows


def create_microstate_review_queue(
    ligand_rows: list[dict[str, Any]],
    output_dir: Path,
) -> list[dict[str, Any]]:
    path = output_dir / "ligand_preparation" / "microstate_review_queue.csv"
    prior = {
        row["ligand_id"]: row
        for row in read_csv(path)
    } if path.exists() else {}
    queue = []
    for row in ligand_rows:
        if not parse_bool(row["state_ambiguity"]):
            continue
        previous = prior.get(row["ligand_id"], {})
        same_selected_state = (
            previous.get("selected_state_smiles", "")
            == row["selected_state_smiles"]
        )
        queue.append({
            "ligand_id": row["ligand_id"],
            "pharmacological_label": row["pharmacological_label"],
            "pH": PH,
            "original_formal_charge": row["original_formal_charge"],
            "selected_formal_charge": row["selected_formal_charge"],
            "selected_state_smiles": row["selected_state_smiles"],
            "alternative_states_json": row["candidate_states_json"],
            "selection_rule": row["protonation_selection_rule"],
            "population_estimated": False,
            "review_status": (
                previous.get("review_status", "PENDING")
                if same_selected_state else "PENDING"
            ),
            "reviewed_state_smiles": (
                previous.get("reviewed_state_smiles", "")
                if same_selected_state else ""
            ),
            "review_note": (
                previous.get("review_note", "") if same_selected_state else ""
            ),
        })
    if queue:
        write_csv(path, queue)
    return queue


def create_mass_docking_preflight(
    ligand_rows: list[dict[str, Any]],
    receptor_preparations: dict[str, dict[str, Any]],
    output_dir: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    for path in (VINA_GPU_EXECUTABLE, VINA_GPU_WRAPPER, VINA_GPU_KERNEL):
        if not path.is_file():
            raise FileNotFoundError(f"Vina-GPU runtime file not found: {path}")
    vina_gpu_version = subprocess.run(
        [str(VINA_GPU_EXECUTABLE), "--version"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    if "Vina-GPU" not in vina_gpu_version:
        raise RuntimeError(f"Unexpected Vina-GPU version output: {vina_gpu_version!r}")
    queue = create_microstate_review_queue(ligand_rows, output_dir)
    reviewed = {
        row["ligand_id"]: row
        for row in queue
        if row["review_status"].strip().upper() == "APPROVED"
        and row["reviewed_state_smiles"].strip()
        == row["selected_state_smiles"].strip()
    }
    manifest = []
    ligand_order = {
        ligand_id: index
        for index, ligand_id in enumerate(sorted(row["ligand_id"] for row in ligand_rows))
    }
    receptor_blockers = {}
    for receptor in RECEPTORS:
        prep = receptor_preparations[receptor]
        low_clearance = [
            row for row in prep["repairs"]
            if row["status"] == "grafted_low_confidence"
        ]
        unresolved = prep["unresolved"]
        receptor_blockers[receptor] = {
            "low_clearance_grafts": [
                {
                    "position": row["uniprot_position"],
                    "clearance_angstrom": row["minimum_new_atom_clearance_angstrom"],
                    "closest_external_atom": row.get("closest_external_atom", ""),
                }
                for row in low_clearance
            ],
            "unresolved_residues": [
                {
                    "position": row.get("uniprot_position", ""),
                    "residue": row.get("source_resname", ""),
                }
                for row in unresolved
            ],
        }
        receptor_path = (
            output_dir / "receptor_preparation" / f"{receptor}_recovered.pdbqt"
        )
        receptor_pdb_path = (
            output_dir / "receptor_preparation" / f"{receptor}_recovered.pdb"
        )
        box_path = BOXES[receptor]
        receptor_ok = not low_clearance and not unresolved
        for ligand in ligand_rows:
            ligand_id = ligand["ligand_id"]
            ligand_path = (
                output_dir / "ligand_preparation" / "pdbqt" / f"{ligand_id}.pdbqt"
            )
            available = (
                parse_bool(ligand["pdbqt_available"])
                and ligand_path.is_file()
            )
            ambiguous = parse_bool(ligand["state_ambiguity"])
            approved = ligand_id in reviewed
            if not available:
                status = "NO DISPONIBLE"
                reason = "Ligando PDBQT preparado ausente"
            elif not receptor_ok:
                status = "REVISIÓN MANUAL"
                reason = "Receptor requiere cierre estructural antes de congelar"
            elif ambiguous and not approved:
                status = "REVISIÓN MANUAL"
                reviewed_row = next(
                    (item for item in queue if item["ligand_id"] == ligand_id),
                    {},
                )
                reason = (
                    "Estado aprobado difiere del protómero preparado; regenerar y validar"
                    if reviewed_row.get("review_status", "").strip().upper() == "APPROVED"
                    else "Estado de microestado pendiente de aprobación dirigida"
                )
            else:
                status = "READY"
                reason = ""
            receptor_index = list(RECEPTORS).index(receptor)
            ligand_index = ligand_order[ligand_id]
            manifest.append(
                {
                    "ligand_id": ligand_id,
                    "receptor": receptor,
                    "status": status,
                    "reason": reason,
                    "ligand_pdbqt": str(ligand_path.relative_to(ROOT)),
                    "ligand_sha256": sha256(ligand_path) if available else "",
                    "receptor_pdbqt": str(receptor_path.relative_to(ROOT)),
                    "receptor_sha256": sha256(receptor_path) if receptor_path.exists() else "",
                    "receptor_pdb": str(receptor_pdb_path.relative_to(ROOT)),
                    "receptor_pdb_sha256": (
                        sha256(receptor_pdb_path) if receptor_pdb_path.exists() else ""
                    ),
                    "box_json": str(box_path.relative_to(ROOT)),
                    "box_sha256": sha256(box_path),
                    "ligand_summary_sha256": sha256(
                        output_dir / "ligand_preparation_summary.csv"
                    ),
                    "recovery_script_sha256": sha256(Path(__file__)),
                    "prolif_script_sha256": sha256(
                        ROOT / "scripts" / "run_prolif_corrected.py"
                    ),
                    "batch_script_sha256": sha256(
                        ROOT / "scripts" / "run_frozen_docking_batch.py"
                    ),
                    "ifp_script_sha256": sha256(
                        ROOT / "scripts" / "run_frozen_ifp.py"
                    ),
                    "vina_executable": str(VINA_GPU_EXECUTABLE),
                    "vina_executable_sha256": sha256(VINA_GPU_EXECUTABLE),
                    "vina_wrapper": str(VINA_GPU_WRAPPER),
                    "vina_wrapper_sha256": sha256(VINA_GPU_WRAPPER),
                    "vina_kernel": str(VINA_GPU_KERNEL),
                    "vina_kernel_sha256": sha256(VINA_GPU_KERNEL),
                    "vina_version": vina_gpu_version,
                    "parameters_json": json.dumps(
                        VINA_GPU_PARAMETERS, sort_keys=True
                    ),
                    "seed": 4513000 + receptor_index * len(ligand_rows) + ligand_index,
                    "pilot_pose_reuse": "NOT_VERIFIED",
                }
            )
    path = output_dir / "mass_docking_preflight.csv"
    write_csv(path, manifest, lineterminator="\n")
    counts = Counter(row["status"] for row in manifest)
    blockers = {
        "microstate_reviews_pending": len(queue) - len(reviewed),
        "receptors": receptor_blockers,
    }
    result = {
        "status": "BLOCKED" if counts["REVISIÓN MANUAL"] else "READY",
        "combinations": len(manifest),
        "status_counts": dict(counts),
        "blockers": blockers,
        "pilot_reuse": (
            "NOT_VERIFIED: pilot runs did not record per-run input hashes; "
            "all combinations will be run in the final batch to avoid stale reuse"
        ),
        "manifest": str(path.relative_to(ROOT)),
    }
    (output_dir / "mass_docking_preflight.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest, result


def write_reports(
    output_dir: Path,
    receptor_preparations: dict[str, dict[str, Any]],
    ligand_rows: list[dict[str, Any]],
    pilot: list[dict[str, Any]],
    docking_rows: list[dict[str, Any]],
    prolif_summary: dict[str, Any],
    decisions: list[dict[str, Any]],
    original_hashes: dict[str, str],
    errors: list[str],
) -> dict[str, Any]:
    repair_rows = [
        row for prep in receptor_preparations.values() for row in prep["repairs"]
    ]
    ligand_decisions = Counter(row["ligand_decision"] for row in ligand_rows)
    docking_decisions = Counter(row["decision"] for row in decisions)
    pilot_ok = sum(row["status"] == "ok" for row in docking_rows)
    prolif_ok = prolif_summary["n_poses"]
    unresolved = [
        (receptor, row["uniprot_position"], row["source_resname"])
        for receptor, prep in receptor_preparations.items()
        for row in prep["unresolved"]
    ]
    low_confidence = [
        (
            row["receptor"],
            row["uniprot_position"],
            row["minimum_new_atom_clearance_angstrom"],
            row.get("closest_external_atom", ""),
        )
        for row in repair_rows
        if row["status"] == "grafted_low_confidence"
    ]
    ambiguous_ligands = [
        row["ligand_id"] for row in ligand_rows
        if parse_bool(row["state_ambiguity"])
    ]
    _, mass_preflight = create_mass_docking_preflight(
        ligand_rows, receptor_preparations, output_dir
    )
    mass_scale_ready = mass_preflight["status"] == "READY"
    conclusion = (
        "APROBADO CON LIMITACIONES DOCUMENTADAS"
        if pilot_ok and prolif_ok and prolif_summary["reproducibility"]["identical"]
        and not unresolved
        else "NO APROBADO"
    )

    receptor_report = [
        "# Preparación localizada de receptores 5-HT2A",
        "",
        "Se preservaron los archivos originales. Solo se transfirieron átomos de cadenas laterales faltantes o de las tres mutaciones de 7WC8 desde receptores experimentales locales, tras ajuste rígido por backbone/CB; no se modeló el resto de la proteína.",
        "",
        "| Receptor | Átomos pesados antes/después | Reparaciones aplicadas | Problemas no reconstruidos | Dictamen |",
        "|---|---:|---:|---|---|",
    ]
    for receptor, prep in receptor_preparations.items():
        applied = sum(row["status"].startswith("grafted") for row in prep["repairs"])
        missing = ", ".join(
            f"{row['source_resname']}{row['source_uniprot_position']}"
            for row in prep["unresolved"]
        ) or "ninguno"
        receptor_report.append(
            f"| {receptor} | {prep['heavy_atoms_before']}/{prep['heavy_atoms_after']} | "
            f"{applied} | {missing} | "
            f"{'preparado para el piloto' if not prep['unresolved'] else 'limitación pendiente'} |"
        )
    receptor_report.extend(
        [
            "",
            "## Residuos injertados y evidencia",
            "",
            "| Receptor | Posición UniProt | Cambio | Plantilla | RMSD backbone/CB (Å) | Clearance mínimo (Å) | Evidencia de selección | Estado |",
            "|---|---:|---|---|---:|---:|---|---|",
        ]
    )
    for row in repair_rows:
        if not row["status"].startswith("grafted"):
            continue
        receptor_report.append(
            f"| {row['receptor']} | {row['uniprot_position']} | "
            f"{row['original_resname']}→{row['prepared_resname']} "
            f"({row['atoms_added']}) | {row['donor_receptor']} | "
            f"{row['local_backbone_cb_rmsd_angstrom']} | "
            f"{row['minimum_new_atom_clearance_angstrom']} | "
            f"{row.get('template_selection_reason', '')} "
            f"{row.get('closest_external_atom', '')} | {row['status']} |"
        )
    receptor_report.extend(
        [
            "",
            "Las sustituciones S162K, M164W y S372N de 7WC8 se trataron como mutaciones ingenierizadas, no como identidad WT. Se generó un modelo de trabajo con SER162, MET164 y SER372 transferidos desde los receptores WT 9LL8/9LLB; se conserva la conformación de backbone de 7WC8. Por ello, este receptor es un modelo de constructo WT-like, no una estructura experimental WT.",
            "",
            "GLU351 de 7WC8 no se reconstruyó: no hay plantilla experimental local completa, su residuo pesado más cercano está a 12.26 Å de las poses piloto y a 4.80 Å fuera del volumen de caja. Se documenta como una limitación periférica; no se inventó su rotámero.",
            "",
            "Para 9LL8 se seleccionó LYS146 de la plantilla experimental 9LLB y ASP218 de 7WC8. La búsqueda combinatoria local mostró que injertar ambos desde 7WC8 deja NZ(LYS146)-OD2(ASP218) a 1.835 Å; la combinación seleccionada aumenta la separación entre esas cadenas laterales a 5.765 Å y conserva clearance mínimo de 2.673 Å frente al resto del receptor. Los ajustes locales backbone/CB fueron 0.170 Å para LYS146 y 0.068 Å para ASP218. Esta corrección cambia el receptor respecto al docking piloto previo de 9LL8; sus 12 poses no se reutilizarán.",
            "",
            "Los campos de RMSD son ajustes locales de backbone/CB y no validan por sí solos la rotámera injertada. Se revisó clearance estérico; el modelo permanece sujeto a incertidumbre conformacional. Los residuos distantes fuera de la envolvente de 6 Å no se tocaron.",
        ]
    )
    (output_dir / "reports" / "receptor_preparation_report.md").write_text(
        "\n".join(receptor_report) + "\n", encoding="utf-8"
    )

    pilot_report = [
        "# Validación piloto de docking y ProLIF 5-HT2A",
        "",
        f"Ligandos elegidos: {len(pilot)} ({sum(item['label'] == 'agonist' for item in pilot)} agonistas, {sum(item['label'] == 'antagonist' for item in pilot)} antagonistas). No se ejecutó el conjunto de 423 ligandos.",
        f"Dockings piloto completados: {pilot_ok}/{len(docking_rows)}. ProLIF: {prolif_ok} modelos de pose. Reproducibilidad de un complejo repetido: {prolif_summary['reproducibility']['identical']}.",
        "",
        "Nota de vigencia: después del piloto se cambió la plantilla de LYS146 de 9LL8 para resolver su contacto estrecho con ASP218. Por ello, los 12 dockings y 84 fingerprints ProLIF previos de 9LL8 son históricos y no validan el receptor actual. No se repetirán como piloto aparte; se recalcularán una sola vez dentro del lote final. Tampoco se reutilizarán poses piloto sin hashes de entrada por ejecución.",
        "",
        "| Receptor | Vina runs OK | Poses ProLIF | Scores Vina mediana (kcal/mol) |",
        "|---|---:|---:|---:|",
    ]
    for receptor in RECEPTORS:
        dr = [row for row in docking_rows if row["receptor"] == receptor]
        scores = [
            float(row["top_affinity_kcal_mol"]) for row in dr
            if row["status"] == "ok" and row["top_affinity_kcal_mol"] != ""
        ]
        pose_count = prolif_summary["receptors"].get(receptor, {}).get("pose_count", 0)
        pilot_report.append(
            f"| {receptor} | {sum(row['status'] == 'ok' for row in dr)}/{len(dr)} | "
            f"{pose_count} | {float(np.median(scores)) if scores else 'n/a'} |"
        )
    pilot_report.extend(
        [
            "",
            "## Frecuencias ProLIF (eventos)",
            "",
            "| Receptor | " + " | ".join(prolif_pipeline.INTERACTIONS) + " |",
            "|---|" + "---:|" * len(prolif_pipeline.INTERACTIONS),
        ]
    )
    for receptor, data in prolif_summary["receptors"].items():
        pilot_report.append(
            f"| {receptor} | "
            + " | ".join(
                str(data["interaction_event_counts"].get(name, 0))
                for name in prolif_pipeline.INTERACTIONS
            )
            + " |"
        )
    pilot_report.extend(
        [
            "",
            "Se conservaron las nueve categorías, incluso cuando la frecuencia es cero. ASP155 se inspeccionó como diagnóstico y no se impuso como filtro. Los scores/poses nuevos no se comparan como equivalentes con los antiguos porque cambiaron el receptor y la protonación de los ligandos; los resultados históricos se conservaron.",
            "",
            f"Parámetros homogéneos del piloto CPU: AutoDock Vina 1.2.7, scoring vina, exhaustiveness={BOX_PARAMETERS['exhaustiveness']}, num_modes={BOX_PARAMETERS['num_modes']}, energy_range={BOX_PARAMETERS['energy_range']}, CPU={BOX_PARAMETERS['cpu']}, seed explícito por combinación. Centro y dimensiones de caja iguales a las configuraciones originales.",
            "",
            f"Controles químicos y del pipeline se registran en `pilot_prolif/pilot_prolif_summary.json`. Errores de preparación/docking: {len(errors)}.",
        ]
    )
    (output_dir / "reports" / "pilot_validation_report.md").write_text(
        "\n".join(pilot_report) + "\n", encoding="utf-8"
    )

    readiness = [
        "# Preparación para continuar el pipeline",
        "",
        f"**Decisión: {conclusion}.**",
        "",
        f"El piloto confirma {pilot_ok}/{len(docking_rows)} dockings, {prolif_ok} poses fingerprintadas y reproducibilidad={prolif_summary['reproducibility']['identical']}. Se mantuvieron separados los tres receptores y las categorías IFP.",
        "",
        f"Ligandos: {len(ligand_rows)} SDF procesados; decisiones de estado: {dict(ligand_decisions)}. El criterio de Dimorphite maximiza protonación de nitrógenos alifáticos no amida y desprotonación de carboxilos a pH 7.4, con estados alternativos registrados. No se estimaron poblaciones. Los casos con ambigüedad se mantienen en revisión manual.",
        "",
        f"Combinaciones de docking: {dict(docking_decisions)}. No se ejecutó redocking masivo. Los resultados nuevos del piloto están en `pilot_docking/`; los resultados originales no se sobrescribieron.",
        "",
        f"Limitaciones estructurales registradas: residuos sin plantilla={unresolved or 'ninguno'}; injertos con clearance pesado mínimo <2 Å={low_confidence or 'ninguno'}. Las conformaciones de cadenas laterales injertadas son modelos de plantilla y deben revisarse; no se interpreta un contacto corto automáticamente como choque sin inspeccionar el par atómico.",
        "",
        f"Estado de escalado masivo: {'LISTO' if mass_scale_ready else 'BLOQUEADO'}. Ligandos con microestado ambiguo pendientes={len(ambiguous_ligands)}. Un estado representativo está preparado para cada uno, pero no se consideró aprobado sin revisión.",
        "",
        f"Preflight actual: {mass_preflight['status']}; {mass_preflight['combinations']} combinaciones incluyendo entradas no disponibles, conteos={mass_preflight['status_counts']}. La cola dirigida está en `ligand_preparation/microstate_review_queue.csv`; solo se aprueba automáticamente el mismo estado ya preparado. Si se selecciona una alternativa, hay que regenerar su PDBQT antes de desbloquear el manifiesto.",
        "",
        "Los 36 dockings del piloto no tienen hashes de entrada por ejecución, así que no se reutilizarán por suposición. El ejecutor final los incluirá una sola vez dentro de las 1.269 combinaciones disponibles; no se harán dockings piloto adicionales.",
        "",
        "La validación ProLIF para el receptor 9LL8 corresponde al modelo anterior al ajuste conjunto LYS146/ASP218; el lote final actualizará docking e IFP en ese receptor.",
        "",
        "El lote final usa Vina-GPU (scoring vina, 1000 lanes, search_depth automático, 9 modos y energy_range=3); no equivale al exhaustiveness=8 del piloto CPU. El manifiesto congela hashes del binario, wrapper y kernel. Tras cerrar la cola, refrescar el preflight y ejecutar `python scripts/run_frozen_docking_batch.py --execute`. Cuando termine, calcular IFP con `python scripts/run_frozen_ifp.py --batch-fingerprint <fingerprint> --execute`.",
        "",
        "El piloto valida el cálculo de fingerprints, no aprueba por sí mismo la geometría de todos los receptores ni resuelve la selección de microestados. No iniciar el lote hasta congelar estas decisiones.",
    ]
    (output_dir / "reports" / "pipeline_readiness_report.md").write_text(
        "\n".join(readiness) + "\n", encoding="utf-8"
    )

    after_hashes = {
        name: sha256(path)
        for name, path in source_paths().items()
        if path.exists()
    }
    changed = [
        name for name, digest in original_hashes.items()
        if after_hashes.get(name) != digest
    ]
    if changed:
        raise RuntimeError(f"Original input hashes changed: {changed}")
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "conclusion": conclusion,
        "tools": {
            "vina": subprocess.run(
                [shutil.which("vina") or "/usr/local/sbin/vina", "--version"],
                capture_output=True, text=True, check=True,
            ).stdout.strip(),
            "openbabel": subprocess.run(
                [shutil.which("obabel") or "/home/sdageri/miniconda3/envs/filter/bin/obabel", "-V"],
                capture_output=True, text=True, check=True,
            ).stderr.strip() or "Open Babel 3.1.0",
            "rdkit": rdBase.rdkitVersion,
            "dimorphite_dl": "2.0.2",
            "prolif": prolif_pipeline.prolif.__version__,
            "pdb2pqr": "3.8.0",
        },
        "parameters": {
            "pH": PH,
            "docking": BOX_PARAMETERS,
            "scoring": "vina",
            "pilot_ligands": len(pilot),
            "mass_redocking_executed": False,
            "full_423_processed": False,
        },
        "receptor_preparations": receptor_preparations,
        "ligand_count": len(ligand_rows),
        "ligand_decisions": dict(ligand_decisions),
        "docking_decisions": dict(docking_decisions),
        "mass_scale_ready": mass_scale_ready,
        "mass_docking_preflight": mass_preflight,
        "mass_scale_blockers": {
            "unresolved_receptor_residues": unresolved,
            "low_clearance_grafts": low_confidence,
            "ambiguous_ligands": ambiguous_ligands,
        },
        "pilot": {
            "ligands": [item["ligand_id"] for item in pilot],
            "docking_runs": len(docking_rows),
            "docking_successes": pilot_ok,
            "prolif_poses": prolif_ok,
            "prolif_reproducibility": prolif_summary["reproducibility"],
            "receptors": prolif_summary["receptors"],
            "validity": {
                receptor: (
                    "historical_only_receptor_updated_after_pilot"
                    if receptor == "9LL8"
                    else "historical_only_input_hashes_not_recorded"
                )
                for receptor in RECEPTORS
            },
        },
        "source_hashes_before_and_after_identical": True,
        "errors": errors,
    }
    (output_dir / "recovery_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return summary


def source_paths() -> dict[str, Path]:
    paths = {
        "ligand_sdf": LIGAND_SDF,
        "selection_csv": SELECTION,
        "receptor_residue_comparison": RESIDUE_MAP,
        "receptor_atom_comparison": ATOM_MAP,
        "ligand_audit": LIGAND_AUDIT,
    }
    for receptor, entry in RECEPTORS.items():
        paths[f"{receptor}_pdb"] = entry["pdb"]
        paths[f"{receptor}_pdbqt"] = PDBQT_INPUTS[receptor]
        paths[f"{receptor}_box"] = BOXES[receptor]
    return paths


def main() -> int:
    if (OUT / "recovery_summary.json").exists():
        raise FileExistsError(
            f"Refusing to overwrite an existing recovery run: {OUT / 'recovery_summary.json'}"
        )
    for path in (LIGAND_SDF, SELECTION, ATOM_MAP, RESIDUE_MAP, LIGAND_AUDIT):
        if not path.exists():
            raise FileNotFoundError(path)
    for receptor, entry in RECEPTORS.items():
        for path in (entry["pdb"], PDBQT_INPUTS[receptor], BOXES[receptor]):
            if not path.exists():
                raise FileNotFoundError(path)
    obabel = shutil.which("obabel") or "/home/sdageri/miniconda3/envs/filter/bin/obabel"
    vina = shutil.which("vina") or "/usr/local/sbin/vina"
    if not Path(obabel).exists() or not Path(vina).exists():
        raise FileNotFoundError("Open Babel and AutoDock Vina are required")
    OUT.mkdir(parents=True, exist_ok=True)
    for subdir in (
        "receptor_preparation", "ligand_preparation", "pilot_docking",
        "pilot_prolif", "logs", "reports",
    ):
        (OUT / subdir).mkdir(parents=True, exist_ok=True)
    original_hashes = {
        name: sha256(path) for name, path in source_paths().items()
    }
    all_residue_rows = read_csv(RESIDUE_MAP)
    all_atom_rows = read_csv(ATOM_MAP)
    receptor_preparations = {}
    for receptor in RECEPTORS:
        summary, _, _ = make_receptor(
            receptor,
            all_residue_rows,
            all_atom_rows,
            all_residue_rows,
            all_atom_rows,
            OUT,
        )
        receptor_preparations[receptor] = summary
    write_csv(
        OUT / "receptor_preparation" / "receptor_repair_manifest.csv",
        [row for item in receptor_preparations.values() for row in item["repairs"]],
    )
    if any(item["unresolved"] for item in receptor_preparations.values()):
        raise RuntimeError(
            "A receptor has an unresolved pocket residue; review receptor_repair_manifest.csv"
        )

    ligand_summary, prepared = prepare_ligands(OUT, obabel)
    pilot = select_pilot(prepared, read_csv(SELECTION), ligand_summary)
    write_csv(
        OUT / "pilot_docking" / "pilot_selection.csv",
        [
            {
                "ligand_id": item["ligand_id"],
                "pharmacological_label": item["label"],
                "scaffold_smiles": item["scaffold"],
                "basic_site_count": item["basic_sites"],
                "hbd": item["hbd"],
                "hba": item["hba"],
                "aromatic_atom_count": item["aromatic_atoms"],
                "heavy_atom_count": item["heavy_atoms"],
                "selected_formal_charge": item["state"]["selected_formal_charge"],
                "state_ambiguity": item["state"]["state_ambiguity"],
            }
            for item in pilot
        ],
    )
    docking_rows = run_pilot_docking(pilot, prepared, OUT, vina)
    errors = [
        f"docking:{row['receptor']}:{row['ligand_id']}:{row['error']}"
        for row in docking_rows if row["status"] != "ok"
    ]
    if errors:
        raise RuntimeError(
            f"Pilot docking had failures: {errors[:5]}; see pilot_docking logs"
        )
    _, prolif_summary = process_pilot_prolif(pilot, prepared, docking_rows, OUT)
    decisions = create_decisions(ligand_summary, receptor_preparations, OUT)
    summary = write_reports(
        OUT, receptor_preparations, ligand_summary, pilot,
        docking_rows, prolif_summary, decisions, original_hashes, errors,
    )
    print(json.dumps({
        "output_dir": str(OUT.relative_to(ROOT)),
        "conclusion": summary["conclusion"],
        "pilot_ligands": len(pilot),
        "docking_runs": len(docking_rows),
        "prolif_poses": summary["pilot"]["prolif_poses"],
        "redocking_decisions": summary["docking_decisions"],
        "mass_redocking_executed": False,
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
