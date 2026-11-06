"""Rebuild and validate the 5-HT2A ProLIF pilot without changing docking inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import prolif
from dimorphite_dl import protonate_smiles
from rdkit import Chem, rdBase
from prolif import Fingerprint, Molecule

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PILOT_REPORT = RESULTS / "prolif" / "pilot" / "pilot_report.json"
AUDIT_REPORT = RESULTS / "prolif" / "audit" / "pilot_report.json"
LIGAND_SDF = ROOT / "structures" / "ligands" / "sdf" / "docking_ligands_selected_3d.sdf"
PDBQT_LIGANDS = ROOT / "structures" / "ligands" / "pdbqt"
INTERACTIONS = (
    "Hydrophobic",
    "HBDonor",
    "HBAcceptor",
    "PiStacking",
    "Anionic",
    "Cationic",
    "CationPi",
    "PiCation",
    "VdWContact",
)
RECEPTORS = {
    "7WC8": {
        "name": "7WC8_chainA_protein",
        "pdb": ROOT / "structures" / "7WC8_chainA_protein.pdb",
        "pdbqt": ROOT / "structures" / "7WC8_chainA_protein.pdbqt",
        "cif": ROOT / "receptors" / "raw" / "7WC8.cif",
    },
    "9LL8": {
        "name": "9LL8_chainR_protein",
        "pdb": ROOT / "structures" / "9LL8_chainR_protein.pdb",
        "pdbqt": ROOT / "structures" / "9LL8_chainR_protein.pdbqt",
        "cif": ROOT / "receptors" / "raw" / "9LL8.cif",
    },
    "9LLB": {
        "name": "9LLB_chainR_protein",
        "pdb": ROOT / "structures" / "9LLB_chainR_protein.pdb",
        "pdbqt": ROOT / "structures" / "9LLB_chainR_protein.pdbqt",
        "cif": ROOT / "receptors" / "raw" / "9LLB.cif",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def atom_residue_key(atom: Chem.Atom) -> tuple[str, int, str, str] | None:
    info = atom.GetPDBResidueInfo()
    if info is None:
        return None
    return (
        info.GetChainId().strip(),
        info.GetResidueNumber(),
        info.GetInsertionCode().strip(),
        info.GetResidueName().strip(),
    )


def atom_key(atom: Chem.Atom) -> tuple[str, int, str, str, str] | None:
    residue = atom_residue_key(atom)
    info = atom.GetPDBResidueInfo()
    if residue is None or info is None:
        return None
    return (*residue, info.GetName().strip())


def pdbqt_atoms(path: Path) -> list[dict[str, Any]]:
    atoms = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_type = line[77:].split()[0] if line[77:].split() else ""
        if not atom_type:
            raise ValueError(f"Missing AutoDock atom type in {path}: {line}")
        elements = {
            "A": "C",
            "C": "C",
            "N": "N",
            "NA": "N",
            "NS": "N",
            "O": "O",
            "OA": "O",
            "S": "S",
            "SA": "S",
            "P": "P",
            "HD": "H",
            "HS": "H",
            "F": "F",
            "Cl": "Cl",
            "CL": "Cl",
            "Br": "Br",
            "BR": "Br",
            "I": "I",
        }
        element = elements.get(atom_type)
        if element is None:
            raise ValueError(f"Unsupported AutoDock atom type {atom_type!r} in {path}")
        try:
            atoms.append(
                {
                    "serial": int(line[6:11]),
                    "atom_name": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21:22].strip(),
                    "resid": int(line[22:26]),
                    "icode": line[26:27].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "charge": float(line[70:76]),
                    "atom_type": atom_type,
                    "element": element,
                }
            )
        except ValueError as exc:
            raise ValueError(f"Malformed PDBQT atom line in {path}: {line}") from exc
    return atoms


def residue_atom_keys_from_pdb(path: Path) -> set[tuple[str, int, str, str, str]]:
    mol = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit could not parse receptor PDB {path}")
    return {
        key
        for atom in mol.GetAtoms()
        if atom.GetSymbol() != "H" and (key := atom_key(atom)) is not None
    }


def audit_docking_pdbqt(receptor_key: str) -> dict[str, Any]:
    entry = RECEPTORS[receptor_key]
    atoms = pdbqt_atoms(entry["pdbqt"])
    heavy = [atom for atom in atoms if atom["element"] != "H"]
    source_heavy = []
    for line in entry["pdb"].read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line[76:78].strip().capitalize()
        if element == "H":
            continue
        source_heavy.append(
            {
                "atom_name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "element": element,
                "chain": line[21].strip(),
                "resid": int(line[22:26]),
                "xyz": tuple(round(float(line[start:start + 8]), 3) for start in (30, 38, 46)),
            }
        )
    source_by_xyz: dict[tuple[float, float, float], list[dict[str, Any]]] = defaultdict(list)
    docking_by_xyz: dict[tuple[float, float, float], list[dict[str, Any]]] = defaultdict(list)
    for atom in source_heavy:
        source_by_xyz[atom["xyz"]].append(atom)
    for atom in heavy:
        atom["xyz"] = tuple(round(atom[axis], 3) for axis in ("x", "y", "z"))
        docking_by_xyz[atom["xyz"]].append(atom)
    common_coordinates = set(source_by_xyz) & set(docking_by_xyz)
    coordinate_matches = sum(
        min(len(source_by_xyz[xyz]), len(docking_by_xyz[xyz]))
        for xyz in common_coordinates
    )
    identity_pairs = []
    for xyz in common_coordinates:
        source_atoms = sorted(
            source_by_xyz[xyz],
            key=lambda atom: (atom["resname"], atom["atom_name"], atom["element"]),
        )
        docking_atoms = sorted(
            docking_by_xyz[xyz],
            key=lambda atom: (atom["resname"], atom["atom_name"], atom["element"]),
        )
        identity_pairs.extend(
            zip(source_atoms, docking_atoms, strict=False)
        )
    atom_name_matches = sum(
        source_atom["atom_name"] == docking_atom["atom_name"]
        for source_atom, docking_atom in identity_pairs
    )
    residue_name_matches = sum(
        source_atom["resname"] == docking_atom["resname"]
        for source_atom, docking_atom in identity_pairs
    )
    element_matches = sum(
        source_atom["element"] == docking_atom["element"]
        for source_atom, docking_atom in identity_pairs
    )
    residue_name_mismatches: dict[tuple[str, int, str, str, int, str], int] = Counter()
    for source_atom, docking_atom in identity_pairs:
        if source_atom["resname"] != docking_atom["resname"]:
            residue_name_mismatches[
                (
                    source_atom["chain"],
                    source_atom["resid"],
                    source_atom["resname"],
                    docking_atom["chain"],
                    docking_atom["resid"],
                    docking_atom["resname"],
                )
            ] += 1
    source_keys = {
        (
            atom["chain"], atom["resid"], atom["resname"], atom["atom_name"]
        )
        for atom in source_heavy
    }
    docking_keys = {
        (
            atom["chain"], atom["resid"], atom["resname"], atom["atom_name"]
        )
        for atom in heavy
    }
    charges = [a["charge"] for a in atoms]
    return {
        "path": str(entry["pdbqt"].relative_to(ROOT)),
        "n_atoms": len(atoms),
        "n_heavy_atoms": len(heavy),
        "n_polar_hydrogens": sum(a["atom_type"] in ("HD", "HS") for a in atoms),
        "atom_type_counts": dict(sorted(Counter(a["atom_type"] for a in atoms).items())),
        "nonzero_partial_charge_atoms": sum(abs(q) > 1e-6 for q in charges),
        "partial_charge_min": min(charges) if charges else None,
        "partial_charge_max": max(charges) if charges else None,
        "heavy_atom_coordinate_matches": coordinate_matches,
        "heavy_atom_coordinate_coverage": coordinate_matches / max(1, len(source_heavy)),
        "atom_name_matches_at_matched_coordinates": atom_name_matches,
        "residue_name_matches_at_matched_coordinates": residue_name_matches,
        "element_matches_at_matched_coordinates": element_matches,
        "source_keys_preserved_verbatim": len(source_keys & docking_keys),
        "residue_chain_or_numbering_was_rewritten": source_keys != docking_keys,
        "residue_name_mismatch_groups": [
            {
                "source_chain": key[0],
                "source_residue_number": key[1],
                "source_residue_name": key[2],
                "pdbqt_chain": key[3],
                "pdbqt_residue_number": key[4],
                "pdbqt_residue_name": key[5],
                "matched_heavy_atoms": count,
            }
            for key, count in sorted(residue_name_mismatches.items())
        ],
    }


def audit_ligand_pdbqt(path: Path) -> dict[str, Any]:
    atoms = pdbqt_atoms(path)
    charges = [atom["charge"] for atom in atoms]
    return {
        "path": str(path.relative_to(ROOT)),
        "n_atoms": len(atoms),
        "n_heavy_atoms": sum(atom["element"] != "H" for atom in atoms),
        "n_polar_hydrogens": sum(atom["atom_type"] in ("HD", "HS") for atom in atoms),
        "atom_type_counts": dict(sorted(Counter(atom["atom_type"] for atom in atoms).items())),
        "nonzero_partial_charge_atoms": sum(abs(charge) > 1e-6 for charge in charges),
        "partial_charge_min": min(charges) if charges else None,
        "partial_charge_max": max(charges) if charges else None,
    }


def residue_mismatch_pose_proximity(
    receptor_key: str,
    mismatch_groups: list[dict[str, Any]],
    baseline: dict[tuple[str, str], dict[str, Any]],
) -> list[dict[str, Any]]:
    source_coords: dict[tuple[str, int, str], list[np.ndarray]] = defaultdict(list)
    for line in RECEPTORS[receptor_key]["pdb"].read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[76:78].strip().capitalize() == "H":
            continue
        source_coords[
            (line[21].strip(), int(line[22:26]), line[17:20].strip())
        ].append(np.asarray([float(line[start:start + 8]) for start in (30, 38, 46)]))
    minimums = {}
    for group in mismatch_groups:
        key = (
            group["source_chain"],
            group["source_residue_number"],
            group["source_residue_name"],
        )
        minimums[key] = {
            "source_chain": key[0],
            "source_residue_number": key[1],
            "source_residue_name": key[2],
            "pdbqt_residue_name": group["pdbqt_residue_name"],
            "minimum_distance_to_pilot_pose_angstrom": math.inf,
            "nearest_ligand_id": None,
            "nearest_pose_model": None,
        }
    if not minimums:
        return []
    for (candidate_receptor, ligand_id), complex_record in baseline.items():
        if candidate_receptor != receptor_key:
            continue
        pose_path = Path(complex_record["pose"])
        if not pose_path.is_absolute():
            pose_path = ROOT / pose_path
        for model in parse_pose_models(pose_path):
            ligand_coords = np.asarray(
                [
                    [atom[axis] for axis in ("x", "y", "z")]
                    for atom in model["atoms"]
                    if atom["element"] != "H"
                ]
            )
            for key, result in minimums.items():
                receptor_coords = np.asarray(source_coords.get(key, []))
                if not len(receptor_coords) or not len(ligand_coords):
                    continue
                distance = float(
                    np.linalg.norm(
                        receptor_coords[:, None, :] - ligand_coords[None, :, :],
                        axis=2,
                    ).min()
                )
                if distance < result["minimum_distance_to_pilot_pose_angstrom"]:
                    result["minimum_distance_to_pilot_pose_angstrom"] = distance
                    result["nearest_ligand_id"] = ligand_id
                    result["nearest_pose_model"] = model["model"]
    return sorted(
        minimums.values(),
        key=lambda item: (
            item["minimum_distance_to_pilot_pose_angstrom"],
            item["source_residue_number"],
        ),
    )


def _pdb2pqr_executable() -> str:
    executable = shutil.which("pdb2pqr")
    if executable is None:
        raise RuntimeError(
            "pdb2pqr is required; run this script in the project filter environment"
        )
    return executable


def prepare_receptor(
    receptor_key: str, pH: float, output_dir: Path
) -> tuple[Molecule, dict[str, Any]]:
    entry = RECEPTORS[receptor_key]
    prep_dir = output_dir / "receptor_preparation"
    prep_dir.mkdir(parents=True, exist_ok=True)
    raw_pdb2pqr = prep_dir / f"{receptor_key}_pdb2pqr_raw.pdb"
    pqr_path = prep_dir / f"{receptor_key}_pH{pH:g}.pqr"
    prepared_path = prep_dir / f"{receptor_key}_prepared.pdb"
    log_path = output_dir / "logs" / f"{receptor_key}_pdb2pqr.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    command = [
        _pdb2pqr_executable(),
        "--ff=AMBER",
        "--titration-state-method=propka",
        f"--with-ph={pH:g}",
        "--keep-chain",
        "--drop-water",
        "--nodebump",
        "--noopt",
        "--pdb-output",
        str(raw_pdb2pqr),
        str(entry["pdb"]),
        str(pqr_path),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode != 0 or not raw_pdb2pqr.exists():
        raise RuntimeError(
            f"PDB2PQR failed for {receptor_key}; see {log_path}: "
            f"{(result.stderr or result.stdout)[-1000:]}"
        )

    source = Chem.MolFromPDBFile(str(entry["pdb"]), sanitize=False, removeHs=False)
    prepared = Chem.MolFromPDBFile(str(raw_pdb2pqr), sanitize=False, removeHs=False)
    if source is None or prepared is None:
        raise RuntimeError(f"RDKit could not parse receptor {receptor_key}")
    source_heavy = {
        key for atom in source.GetAtoms()
        if atom.GetSymbol() != "H" and (key := atom_key(atom)) is not None
    }
    source_coords = {
        atom_key(atom): np.asarray(source.GetConformer().GetAtomPosition(atom.GetIdx()))
        for atom in source.GetAtoms()
        if atom.GetSymbol() != "H" and atom_key(atom) is not None
    }
    conf = prepared.GetConformer()
    retained_heavy = {
        atom.GetIdx()
        for atom in prepared.GetAtoms()
        if atom.GetSymbol() != "H" and atom_key(atom) in source_heavy
    }
    if len(retained_heavy) != len(source_heavy):
        raise RuntimeError(
            f"PDB2PQR did not preserve the source heavy-atom inventory for {receptor_key}"
        )
    retained = set(retained_heavy)
    for atom in prepared.GetAtoms():
        if atom.GetSymbol() != "H":
            continue
        candidates = []
        for heavy_index in retained_heavy:
            heavy_atom = prepared.GetAtomWithIdx(heavy_index)
            if atom_residue_key(atom) != atom_residue_key(heavy_atom):
                continue
            distance = (
                conf.GetAtomPosition(atom.GetIdx())
                - conf.GetAtomPosition(heavy_index)
            ).Length()
            if distance <= 1.3:
                candidates.append((distance, heavy_index))
        if candidates:
            retained.add(atom.GetIdx())

    rw = Chem.RWMol(prepared)
    for index in sorted(set(range(prepared.GetNumAtoms())) - retained, reverse=True):
        rw.RemoveAtom(index)
    mol = rw.GetMol()

    removed_cross_residue_bonds = []
    conf = mol.GetConformer()
    rw = Chem.RWMol(mol)
    for bond in list(mol.GetBonds()):
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        residue1 = atom_residue_key(atom1)
        residue2 = atom_residue_key(atom2)
        if residue1 == residue2:
            continue
        info1 = atom1.GetPDBResidueInfo()
        info2 = atom2.GetPDBResidueInfo()
        name1 = info1.GetName().strip() if info1 else ""
        name2 = info2.GetName().strip() if info2 else ""
        resname1 = info1.GetResidueName().strip() if info1 else ""
        resname2 = info2.GetResidueName().strip() if info2 else ""
        distance = (
            conf.GetAtomPosition(atom1.GetIdx())
            - conf.GetAtomPosition(atom2.GetIdx())
        ).Length()
        peptide = (
            residue1 is not None
            and residue2 is not None
            and residue1[0] == residue2[0]
            and distance <= 1.8
            and (
                (name1 == "C" and name2 == "N")
                or (name1 == "N" and name2 == "C")
            )
        )
        disulfide = (
            resname1 == resname2 == "CYS"
            and name1 == name2 == "SG"
            and 1.8 <= distance <= 2.3
        )
        if not (peptide or disulfide):
            removed_cross_residue_bonds.append(
                {
                    "atom1": atom_key(atom1),
                    "atom2": atom_key(atom2),
                    "distance": round(distance, 3),
                }
            )
            rw.RemoveBond(atom1.GetIdx(), atom2.GetIdx())
    mol = rw.GetMol()

    rw = Chem.RWMol(mol)
    added_hydrogen_bonds = []
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "H":
            continue
        candidates = []
        for heavy_atom in mol.GetAtoms():
            if heavy_atom.GetSymbol() == "H":
                continue
            if atom_residue_key(atom) != atom_residue_key(heavy_atom):
                continue
            distance = (
                conf.GetAtomPosition(atom.GetIdx())
                - conf.GetAtomPosition(heavy_atom.GetIdx())
            ).Length()
            if distance <= 1.3:
                candidates.append((distance, heavy_atom.GetIdx()))
        if not candidates:
            continue
        _, parent = min(candidates)
        if mol.GetBondBetweenAtoms(atom.GetIdx(), parent) is None:
            rw.AddBond(atom.GetIdx(), parent, Chem.BondType.SINGLE)
            added_hydrogen_bonds.append(
                {"hydrogen": atom_key(atom), "parent": atom_key(mol.GetAtomWithIdx(parent))}
            )
    mol = rw.GetMol()

    carboxylate_adjustments = []
    acid_groups: dict[tuple[str, int, str, str], dict[str, int]] = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        residue = atom_residue_key(atom)
        if info is None or residue is None:
            continue
        resname = info.GetResidueName().strip()
        atom_name = info.GetName().strip()
        if (
            (resname == "ASP" and atom_name in ("CG", "OD1", "OD2"))
            or (resname == "GLU" and atom_name in ("CD", "OE1", "OE2"))
        ):
            acid_groups.setdefault(residue, {})[atom_name] = atom.GetIdx()
    for residue, atoms in acid_groups.items():
        resname = residue[3]
        carbon_name, oxygen_names = (
            ("CG", ("OD1", "OD2"))
            if resname == "ASP"
            else ("CD", ("OE1", "OE2"))
        )
        carbon = atoms.get(carbon_name)
        oxygens = [atoms.get(name) for name in oxygen_names]
        if carbon is None or any(index is None for index in oxygens):
            continue
        protonated = [
            index
            for index in oxygens
            if any(neighbor.GetSymbol() == "H" for neighbor in mol.GetAtomWithIdx(index).GetNeighbors())
        ]
        negatively_charged = [
            index for index in oxygens
            if mol.GetAtomWithIdx(index).GetFormalCharge() < 0
        ]
        single_bond_oxygen = (
            protonated[0]
            if protonated
            else (negatively_charged[0] if negatively_charged else oxygens[0])
        )
        for index in oxygens:
            bond = mol.GetBondBetweenAtoms(carbon, index)
            if bond is None:
                raise RuntimeError(
                    f"Missing {resname} carboxyl bond at {residue}; refusing to infer it"
                )
            old_type = str(bond.GetBondType())
            desired_type = (
                Chem.BondType.SINGLE
                if index == single_bond_oxygen
                else Chem.BondType.DOUBLE
            )
            bond.SetBondType(desired_type)
            mol.GetAtomWithIdx(index).SetFormalCharge(
                0 if protonated else (-1 if index == single_bond_oxygen else 0)
            )
            if old_type != str(desired_type):
                carboxylate_adjustments.append(
                    {
                        "residue": residue,
                        "atom": atom_key(mol.GetAtomWithIdx(index)),
                        "old_bond_order": old_type,
                        "new_bond_order": str(desired_type),
                        "protonation": "neutral_acid" if protonated else "carboxylate",
                    }
                )

    problems = Chem.DetectChemistryProblems(mol)
    if problems:
        details = []
        for problem in problems:
            index = problem.GetAtomIdx()
            atom = mol.GetAtomWithIdx(index)
            details.append(
                f"{problem.Message()} at {atom_key(atom)}"
            )
        raise RuntimeError(
            f"Receptor chemistry remains invalid for {receptor_key}: "
            + "; ".join(details[:10])
        )
    Chem.SanitizeMol(mol)

    prepared_coords = {
        atom_key(atom): np.asarray(mol.GetConformer().GetAtomPosition(atom.GetIdx()))
        for atom in mol.GetAtoms()
        if atom.GetSymbol() != "H" and atom_key(atom) is not None
    }
    if prepared_coords.keys() != source_coords.keys():
        raise RuntimeError(f"Heavy-atom identity set changed for {receptor_key}")
    max_heavy_coordinate_delta = max(
        float(np.linalg.norm(prepared_coords[key] - source_coords[key]))
        for key in source_coords
    )
    if max_heavy_coordinate_delta > 0.001:
        raise RuntimeError(
            f"Heavy-atom coordinates moved for {receptor_key}: "
            f"{max_heavy_coordinate_delta:.4f} A"
        )

    Chem.MolToPDBFile(mol, str(prepared_path))
    protein = Molecule.from_rdkit(mol)
    receptor_report = {
        "receptor_key": receptor_key,
        "source_pdb": str(entry["pdb"].relative_to(ROOT)),
        "source_cif": str(entry["cif"].relative_to(ROOT)),
        "source_docking_pdbqt": str(entry["pdbqt"].relative_to(ROOT)),
        "prepared_pdb": str(prepared_path.relative_to(ROOT)),
        "pH": pH,
        "protonation_method": "PDB2PQR 3.8.0 with PROPKA; AMBER force field; nodebump/noopt",
        "water_policy": "drop-water",
        "nonstandard_residue_policy": "PDB2PQR-supported residues only; unsupported residues fail explicitly",
        "source_heavy_atoms": len(source_heavy),
        "prepared_heavy_atoms": sum(atom.GetSymbol() != "H" for atom in mol.GetAtoms()),
        "explicit_hydrogens": sum(atom.GetSymbol() == "H" for atom in mol.GetAtoms()),
        "pdb2pqr_added_heavy_atoms_discarded": len(
            [
                atom for atom in prepared.GetAtoms()
                if atom.GetSymbol() != "H" and atom_key(atom) not in source_heavy
            ]
        ),
        "removed_nonlocal_inter_residue_bonds": removed_cross_residue_bonds,
        "added_same_residue_hydrogen_bonds": added_hydrogen_bonds,
        "carboxylate_bond_order_adjustments": carboxylate_adjustments,
        "diagnostic_residues": inspect_diagnostic_residues(mol),
        "max_heavy_atom_coordinate_delta_angstrom": max_heavy_coordinate_delta,
        "sanitization_problems_after_corrections": 0,
        "pdb2pqr_log": str(log_path.relative_to(ROOT)),
    }
    (prep_dir / f"{receptor_key}_preparation.json").write_text(
        json.dumps(receptor_report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return protein, receptor_report


def read_sdf_records(path: Path) -> dict[str, Chem.Mol]:
    records = {}
    for mol in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True):
        if mol is None or not mol.HasProp("ligand_id"):
            continue
        ligand_id = mol.GetProp("ligand_id")
        if ligand_id in records:
            raise ValueError(f"Duplicate ligand_id {ligand_id} in {path}")
        records[ligand_id] = mol
    return records


def heavy_graph(mol: Chem.Mol) -> tuple[nx.Graph, list[int]]:
    heavy_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() != "H"]
    graph = nx.Graph()
    heavy_set = set(heavy_indices)
    for index in heavy_indices:
        atom = mol.GetAtomWithIdx(index)
        degree = sum(neighbor.GetIdx() in heavy_set for neighbor in atom.GetNeighbors())
        graph.add_node(
            index,
            atomic_number=atom.GetAtomicNum(),
            degree=degree,
            aromatic=atom.GetIsAromatic(),
        )
    for bond in mol.GetBonds():
        first = bond.GetBeginAtomIdx()
        second = bond.GetEndAtomIdx()
        if first in heavy_set and second in heavy_set:
            graph.add_edge(first, second)
    return graph, heavy_indices


def _kabsch_rmsd(
    first: np.ndarray, second: np.ndarray
) -> float:
    first_centered = first - first.mean(axis=0)
    second_centered = second - second.mean(axis=0)
    u, _, vt = np.linalg.svd(first_centered.T @ second_centered)
    determinant = np.linalg.det(u @ vt)
    rotation = u @ np.diag([1.0, 1.0, determinant]) @ vt
    delta = first_centered @ rotation - second_centered
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def ligand_atom_mapping(
    reference: Chem.Mol, prepared_pdbqt: Path, obabel: str, scratch: Path
) -> tuple[Chem.Mol, dict[int, int], dict[str, Any]]:
    converted_sdf = scratch / f"{prepared_pdbqt.stem}_recovered.sdf"
    process = subprocess.run(
        [obabel, "-ipdbqt", str(prepared_pdbqt), "-osdf", "-O", str(converted_sdf)],
        capture_output=True,
        text=True,
    )
    if process.returncode != 0 or not converted_sdf.exists():
        raise RuntimeError(
            f"Open Babel could not inspect {prepared_pdbqt}: "
            f"{process.stderr[-500:]}"
        )
    converted = Chem.SDMolSupplier(str(converted_sdf), removeHs=False, sanitize=True)[0]
    if converted is None:
        raise RuntimeError(f"Open Babel output is unreadable: {converted_sdf}")

    reference_graph, reference_indices = heavy_graph(reference)
    converted_graph, converted_indices = heavy_graph(converted)
    if len(reference_indices) != len(converted_indices):
        raise ValueError(
            f"Heavy-atom count mismatch: SDF={len(reference_indices)}, "
            f"prepared PDBQT={len(converted_indices)}"
        )
    node_match = lambda left, right: all(
        left[name] == right[name]
        for name in ("atomic_number", "degree", "aromatic")
    )
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        reference_graph, converted_graph, node_match=node_match
    )
    reference_xyz = reference.GetConformer().GetPositions()
    converted_xyz = converted.GetConformer().GetPositions()
    ranked = []
    for mapping in matcher.isomorphisms_iter():
        ref_order = sorted(mapping)
        rmsd = _kabsch_rmsd(
            np.asarray([reference_xyz[index] for index in ref_order]),
            np.asarray([converted_xyz[mapping[index]] for index in ref_order]),
        )
        ranked.append((rmsd, mapping))
        if len(ranked) > 10000:
            raise ValueError("Too many symmetric atom mappings to establish identity")
    ranked.sort(key=lambda item: (item[0], tuple(sorted(item[1].items()))))
    if not ranked:
        raise ValueError("No topology-preserving atom mapping from SDF to prepared PDBQT")
    best_rmsd = ranked[0][0]
    second_rmsd = ranked[1][0] if len(ranked) > 1 else None
    if best_rmsd > 0.05:
        raise ValueError(f"Best reference conformer mapping RMSD is {best_rmsd:.3f} A")
    if second_rmsd is not None and second_rmsd - best_rmsd < 0.05:
        raise ValueError(
            f"Atom mapping is geometrically ambiguous: best={best_rmsd:.3f} A, "
            f"next={second_rmsd:.3f} A"
        )
    mapping = ranked[0][1]

    prepared_atoms = pdbqt_atoms(prepared_pdbqt)
    pdbqt_heavy = [atom for atom in prepared_atoms if atom["element"] != "H"]
    converted_heavy = [
        index for index, atom in enumerate(converted.GetAtoms())
        if atom.GetSymbol() != "H"
    ]
    if len(pdbqt_heavy) != len(converted_heavy):
        raise ValueError("PDBQT heavy-atom list differs from Open Babel reconstruction")
    for record, index in zip(pdbqt_heavy, converted_heavy, strict=True):
        atom = converted.GetAtomWithIdx(index)
        if record["element"] != atom.GetSymbol():
            raise ValueError(
                f"Element mismatch at PDBQT serial {record['serial']}: "
                f"{record['element']} != {atom.GetSymbol()}"
            )
        xyz = np.asarray([record["x"], record["y"], record["z"]])
        if float(np.linalg.norm(xyz - converted_xyz[index])) > 0.002:
            raise ValueError("PDBQT/Open Babel atom order is not coordinate-identical")

    return converted, mapping, {
        "reference_heavy_atoms": len(reference_indices),
        "pdbqt_heavy_atoms": len(pdbqt_heavy),
        "graph_isomorphisms": len(ranked),
        "best_alignment_rmsd_angstrom": best_rmsd,
        "next_alignment_rmsd_angstrom": second_rmsd,
        "mapping_status": "unique_by_graph_and_conformer_geometry",
        "sdf_to_pdbqt_atom_map": [
            {
                "sdf_atom_index": reference_index,
                "pdbqt_serial": pdbqt_heavy[
                    converted_heavy.index(converted_index)
                ]["serial"],
                "element": converted.GetAtomWithIdx(converted_index).GetSymbol(),
            }
            for reference_index, converted_index in sorted(mapping.items())
        ],
    }


def _is_carbonyl_carbon(atom: Chem.Atom) -> bool:
    if atom.GetSymbol() != "C":
        return False
    return any(
        bond.GetBondType() == Chem.BondType.DOUBLE
        and bond.GetOtherAtom(atom).GetSymbol() in ("O", "S")
        for bond in atom.GetBonds()
    )


def inspect_diagnostic_residues(mol: Chem.Mol) -> dict[str, Any]:
    expected_atoms = {
        "ASP155": ("ASP", 155, {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"}),
        "PHE339": ("PHE", 339, {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
        "PHE340": ("PHE", 340, {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
        "ASN343": ("ASN", 343, {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"}),
        "TRP336": ("TRP", 336, {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}),
    }
    residues: dict[tuple[str, int, str, str], list[Chem.Atom]] = defaultdict(list)
    for atom in mol.GetAtoms():
        residue = atom_residue_key(atom)
        if residue is not None:
            residues[residue].append(atom)
    result = {}
    for label, (resname, number, expected) in expected_atoms.items():
        matches = [
            (key, atoms)
            for key, atoms in residues.items()
            if key[1] == number and key[3] == resname
        ]
        if len(matches) > 1:
            raise RuntimeError(f"Ambiguous diagnostic residue {label}: {matches}")
        if not matches:
            result[label] = {"present": False}
            continue
        key, atoms = matches[0]
        heavy_atoms = [atom for atom in atoms if atom.GetSymbol() != "H"]
        names = {atom.GetPDBResidueInfo().GetName().strip() for atom in heavy_atoms}
        result[label] = {
            "present": True,
            "chain": key[0],
            "residue_number": key[1],
            "residue_name": key[3],
            "heavy_atom_names": sorted(names),
            "missing_expected_heavy_atoms": sorted(expected - names),
            "aromatic_heavy_atoms": sorted(
                atom.GetPDBResidueInfo().GetName().strip()
                for atom in heavy_atoms
                if atom.GetIsAromatic()
            ),
            "explicit_hydrogen_count": sum(
                atom.GetSymbol() == "H" for atom in atoms
            ),
            "carboxyl_oxygen_formal_charges": {
                atom.GetPDBResidueInfo().GetName().strip(): atom.GetFormalCharge()
                for atom in heavy_atoms
                if label == "ASP155"
                and atom.GetPDBResidueInfo().GetName().strip() in {"OD1", "OD2"}
            },
        }
    return result


def _basic_nitrogens(mol: Chem.Mol) -> set[int]:
    basic = set()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N" or atom.GetIsAromatic():
            continue
        heavy_neighbors = [
            neighbor for neighbor in atom.GetNeighbors()
            if neighbor.GetSymbol() != "H"
        ]
        if not any(neighbor.GetSymbol() == "C" for neighbor in heavy_neighbors):
            continue
        if any(_is_carbonyl_carbon(neighbor) for neighbor in heavy_neighbors):
            continue
        if any(
            bond.GetBondType() == Chem.BondType.DOUBLE
            for bond in atom.GetBonds()
        ):
            continue
        basic.add(atom.GetIdx())
    return basic


def _acid_oxygens(mol: Chem.Mol) -> set[int]:
    acidic = set()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O" or atom.GetIsAromatic():
            continue
        for neighbor in atom.GetNeighbors():
            if _is_carbonyl_carbon(neighbor):
                acidic.add(atom.GetIdx())
    return acidic


def select_protomer(
    reference: Chem.Mol, pH: float
) -> tuple[Chem.Mol, dict[int, int], dict[str, Any]]:
    smiles = Chem.MolToSmiles(
        Chem.RemoveHs(Chem.Mol(reference)), canonical=True, isomericSmiles=True
    )
    candidates = protonate_smiles(
        smiles,
        ph_min=pH,
        ph_max=pH,
        precision=0.25,
        max_variants=128,
        validate_output=True,
    )
    if not candidates:
        raise ValueError(f"Dimorphite returned no protonation state at pH {pH:g}")
    reference_heavy = Chem.RemoveHs(Chem.Mol(reference))
    ref_graph, ref_indices = heavy_graph(reference_heavy)
    basic_atoms = _basic_nitrogens(reference_heavy)
    acid_atoms = _acid_oxygens(reference_heavy)
    states = []
    for smiles_candidate in sorted(set(candidates)):
        candidate = Chem.MolFromSmiles(smiles_candidate)
        if candidate is None:
            raise ValueError(f"Dimorphite returned invalid SMILES: {smiles_candidate}")
        cand_graph, cand_indices = heavy_graph(candidate)
        if len(cand_indices) != len(ref_indices):
            raise ValueError("Dimorphite changed the heavy-atom count")
        matcher = nx.algorithms.isomorphism.GraphMatcher(
            ref_graph,
            cand_graph,
            node_match=lambda left, right: all(
                left[name] == right[name]
                for name in ("atomic_number", "degree", "aromatic")
            ),
        )
        mappings = list(matcher.isomorphisms_iter())
        if not mappings:
            raise ValueError("Dimorphite changed heavy-atom connectivity")
        ranked_maps = []
        for mapping in mappings:
            score = sum(
                candidate.GetAtomWithIdx(mapping[index]).GetFormalCharge() > 0
                for index in basic_atoms
            ) + sum(
                candidate.GetAtomWithIdx(mapping[index]).GetFormalCharge() < 0
                for index in acid_atoms
            )
            ranked_maps.append((score, tuple(sorted(mapping.items())), mapping))
        ranked_maps.sort(key=lambda item: (-item[0], item[1]))
        score, _, mapping = ranked_maps[0]
        states.append(
            {
                "smiles": smiles_candidate,
                "mol": candidate,
                "mapping": mapping,
                "score": score,
                "formal_charge": Chem.GetFormalCharge(candidate),
                "mapping_count": len(mappings),
            }
        )
    states.sort(
        key=lambda state: (
            -state["score"],
            -state["formal_charge"],
            state["smiles"],
        )
    )
    chosen = states[0]
    same_score = [state for state in states if state["score"] == chosen["score"]]
    unique_states = {state["smiles"] for state in same_score}
    ambiguous = len(unique_states) > 1
    return chosen["mol"], chosen["mapping"], {
        "input_smiles": smiles,
        "pH": pH,
        "candidate_states": len(states),
        "chosen_smiles": chosen["smiles"],
        "formal_charge": chosen["formal_charge"],
        "basic_sites_protonated_score": chosen["score"],
        "basic_site_count": len(basic_atoms),
        "carboxyl_oxygen_count": len(acid_atoms),
        "state_ambiguity": ambiguous,
        "candidate_smiles": [state["smiles"] for state in states],
        "selection_policy": (
            "maximize protonation of non-amide aliphatic N sites and "
            "deprotonation of carboxyl groups; tie-break by net charge, "
            "then canonical isomeric SMILES"
        ),
        "warning": (
            "Dimorphite enumerates plausible states but does not provide "
            "populations; selected state is a deterministic model, not a "
            "measured dominant microstate."
            if ambiguous
            else ""
        ),
    }


def parse_pose_models(path: Path) -> list[dict[str, Any]]:
    models = []
    current: list[str] | None = None
    current_number = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL "):
            current_number = int(line.split()[1])
            current = []
        elif line.startswith("ENDMDL"):
            if current is not None:
                models.append(
                    {"model": current_number, "atoms": parse_pdbqt_atom_lines(current)}
                )
            current = None
        elif current is not None:
            current.append(line)
    if not models:
        lines = path.read_text(encoding="utf-8").splitlines()
        models = [{"model": 1, "atoms": parse_pdbqt_atom_lines(lines)}]
    return models


def parse_pdbqt_atom_lines(lines: list[str]) -> list[dict[str, Any]]:
    atoms = []
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_type = line[77:].split()[0] if line[77:].split() else ""
        element = {
            "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
            "O": "O", "OA": "O", "S": "S", "SA": "S", "P": "P",
            "HD": "H", "HS": "H", "F": "F", "Cl": "Cl", "CL": "Cl",
            "Br": "Br", "BR": "Br", "I": "I",
        }.get(atom_type)
        if element is None:
            raise ValueError(f"Unsupported pose atom type {atom_type!r}")
        atoms.append(
            {
                "serial": int(line[6:11]),
                "atom_name": line[12:16].strip(),
                "atom_type": atom_type,
                "element": element,
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            }
        )
    return atoms


def reconstruct_pose_ligand(
    reference: Chem.Mol,
    protomer: Chem.Mol,
    protomer_mapping: dict[int, int],
    prepared_atom_mapping: dict[int, int],
    prepared_heavy_indices: list[int],
    pose_atoms: list[dict[str, Any]],
) -> tuple[Chem.Mol, float]:
    pose_heavy = [atom for atom in pose_atoms if atom["element"] != "H"]
    pose_hydrogens = [atom for atom in pose_atoms if atom["element"] == "H"]
    if len(pose_heavy) != len(prepared_heavy_indices):
        raise ValueError(
            f"Pose has {len(pose_heavy)} heavy atoms; template has "
            f"{len(prepared_heavy_indices)}"
        )
    for pose_atom, prepared_index in zip(
        pose_heavy, prepared_heavy_indices, strict=True
    ):
        prepared_atom = prepared_index
        if pose_atom["element"] != prepared_atom_mapping["elements"][prepared_atom]:
            raise ValueError("Pose/PDBQT heavy atom element ordering differs")
    # prepared_atom_mapping carries the validated reference-SDF -> PDBQT index map.
    prepared_order = {
        atom_index: order for order, atom_index in enumerate(prepared_heavy_indices)
    }
    ref_to_pose = {
        ref_index: prepared_order[prepared_index]
        for ref_index, prepared_index in prepared_atom_mapping["ref_to_prepared"].items()
    }
    conformer = Chem.Conformer(protomer.GetNumAtoms())
    for reference_index, protomer_index in protomer_mapping.items():
        pose_atom = pose_heavy[ref_to_pose[reference_index]]
        conformer.SetAtomPosition(
            protomer_index, (pose_atom["x"], pose_atom["y"], pose_atom["z"])
        )
    ligand = Chem.Mol(protomer)
    ligand.RemoveAllConformers()
    ligand.AddConformer(conformer, assignId=True)
    ligand = Chem.AddHs(ligand, addCoords=True)

    # Retain docking hydrogens where the PDBQT pose contains a polar H.
    heavy_lookup = []
    for reference_index, protomer_index in protomer_mapping.items():
        point = np.asarray(
            ligand.GetConformer().GetAtomPosition(protomer_index), dtype=float
        )
        heavy_lookup.append((reference_index, protomer_index, point))
    assigned_h = set()
    conf = ligand.GetConformer()
    for hydrogen in pose_hydrogens:
        point = np.asarray([hydrogen["x"], hydrogen["y"], hydrogen["z"]])
        parent_candidates = sorted(
            (
                float(np.linalg.norm(point - heavy_point)),
                reference_index,
                protomer_index,
            )
            for reference_index, protomer_index, heavy_point in heavy_lookup
        )
        if not parent_candidates or parent_candidates[0][0] > 1.3:
            continue
        _, _, parent_index = parent_candidates[0]
        ligand_hydrogens = [
            neighbor.GetIdx()
            for neighbor in ligand.GetAtomWithIdx(parent_index).GetNeighbors()
            if neighbor.GetSymbol() == "H" and neighbor.GetIdx() not in assigned_h
        ]
        if ligand_hydrogens:
            hydrogen_index = ligand_hydrogens[0]
            conf.SetAtomPosition(
                hydrogen_index, (hydrogen["x"], hydrogen["y"], hydrogen["z"])
            )
            assigned_h.add(hydrogen_index)

    max_coordinate_delta = 0.0
    prepared_order = {
        atom_index: order
        for order, atom_index in enumerate(prepared_heavy_indices)
    }
    for reference_index, protomer_index in protomer_mapping.items():
        prepared_index = prepared_atom_mapping["ref_to_prepared"][reference_index]
        pose_index = prepared_order[prepared_index]
        expected = np.asarray(
            [pose_heavy[pose_index][axis] for axis in ("x", "y", "z")]
        )
        observed = np.asarray(conf.GetAtomPosition(protomer_index))
        max_coordinate_delta = max(
            max_coordinate_delta, float(np.linalg.norm(expected - observed))
        )
    if max_coordinate_delta > 0.001:
        raise ValueError(
            f"Heavy-atom pose coordinates changed by {max_coordinate_delta:.4f} A"
        )
    Chem.SanitizeMol(ligand)
    return ligand, max_coordinate_delta


def baseline_models(report_path: Path) -> dict[tuple[str, str], dict[str, Any]]:
    data = json.loads(report_path.read_text(encoding="utf-8"))
    result = {}
    for complex_record in data.get("complexes", []):
        key = (complex_record["receptor_key"], complex_record["ligand_id"])
        if key in result:
            raise ValueError(f"Duplicate complex in baseline report: {key}")
        result[key] = complex_record
    return result


def baseline_counts(complexes: dict[tuple[str, str], dict[str, Any]]) -> dict[str, int]:
    counts = Counter()
    for complex_record in complexes.values():
        for model in complex_record.get("models", []):
            for event in model.get("interactions", []):
                counts[event.get("interaction", "")] += 1
    return {name: counts.get(name, 0) for name in INTERACTIONS}


def interaction_rows(ifp: dict) -> tuple[list[dict[str, Any]], set[str]]:
    rows = []
    binary_features = set()
    for (ligand_residue, protein_residue), terms in ifp.items():
        protein_id = str(protein_residue)
        for interaction, metadata_list in terms.items():
            for metadata in metadata_list:
                rows.append(
                    {
                        "ligand_residue": str(ligand_residue),
                        "protein_residue": protein_id,
                        "interaction": str(interaction),
                        "distance": json_safe(metadata.get("distance")),
                        "indices": json_safe(metadata.get("indices")),
                    }
                )
                binary_features.add(
                    f"{ligand_residue}|{protein_id}|{interaction}"
                )
    rows.sort(
        key=lambda row: (
            row["protein_residue"],
            row["interaction"],
            json.dumps(row["indices"], sort_keys=True),
        )
    )
    return rows, binary_features


def feature_schema(protein: Molecule, ligand: Molecule) -> list[str]:
    protein_residues = sorted(str(key) for key in protein.residues)
    ligand_residues = sorted(str(key) for key in ligand.residues)
    return [
        f"{ligand_residue}|{protein_residue}|{interaction}"
        for ligand_residue in ligand_residues
        for protein_residue in protein_residues
        for interaction in INTERACTIONS
    ]


def process_one_pose(
    receptor_key: str,
    complex_record: dict[str, Any],
    model_number: int,
    pose_atoms: list[dict[str, Any]],
    reference: Chem.Mol,
    protomer: Chem.Mol,
    protomer_mapping: dict[int, int],
    prepared_atom_mapping: dict[str, Any],
    protein: Molecule,
    output_dir: Path,
) -> dict[str, Any]:
    ligand, coordinate_delta = reconstruct_pose_ligand(
        reference,
        protomer,
        protomer_mapping,
        prepared_atom_mapping,
        prepared_atom_mapping["prepared_heavy_indices"],
        pose_atoms,
    )
    ligand_molecule = Molecule.from_rdkit(ligand)
    fingerprint = Fingerprint(
        interactions=list(INTERACTIONS),
        count=False,
        implicit_hydrogens=False,
    )
    ifp = fingerprint.generate(ligand_molecule, protein, metadata=True)
    events, binary_features = interaction_rows(ifp)
    counts = {name: 0 for name in INTERACTIONS}
    for event in events:
        counts[event["interaction"]] += 1
    binary_counts = {
        name: int(any(event["interaction"] == name for event in events))
        for name in INTERACTIONS
    }
    identifier = (
        f"{receptor_key}_{complex_record['ligand_id']}_model_{model_number:02d}"
    )
    ligand_dir = output_dir / "ligand_reconstruction" / receptor_key
    complex_dir = output_dir / "complexes" / receptor_key
    fingerprint_dir = output_dir / "fingerprints" / receptor_key
    for directory in (ligand_dir, complex_dir, fingerprint_dir):
        directory.mkdir(parents=True, exist_ok=True)
    ligand_sdf = ligand_dir / f"{identifier}.sdf"
    writer = Chem.SDWriter(str(ligand_sdf))
    ligand.SetProp("ligand_id", complex_record["ligand_id"])
    ligand.SetProp("receptor_key", receptor_key)
    ligand.SetProp("pose_model", str(model_number))
    ligand.SetProp(
        "pharmacological_label",
        complex_record.get("pharmacological_label", ""),
    )
    writer.write(ligand)
    writer.close()

    residue_keys = feature_schema(protein, ligand_molecule)
    bits = np.asarray([int(feature in binary_features) for feature in residue_keys], dtype=np.uint8)
    fingerprint_path = fingerprint_dir / f"{identifier}.npz"
    np.savez_compressed(
        fingerprint_path,
        bits=bits,
        feature_names=np.asarray(residue_keys, dtype=str),
        interaction_types=np.asarray(INTERACTIONS, dtype=str),
    )
    payload = {
        "receptor_key": receptor_key,
        "ligand_id": complex_record["ligand_id"],
        "pharmacological_label": complex_record.get("pharmacological_label", ""),
        "pose_model": model_number,
        "pose_path": complex_record["pose"],
        "pose_sha256": sha256(Path(complex_record["pose"])),
        "receptor_prepared_pdb": (
            f"receptor_preparation/{receptor_key}_prepared.pdb"
        ),
        "status": "processed",
        "n_interaction_events": len(events),
        "interaction_event_counts": counts,
        "binary_category_presence": binary_counts,
        "max_heavy_atom_coordinate_delta_angstrom": coordinate_delta,
        "ligand_sdf": str(ligand_sdf.relative_to(ROOT)),
        "fingerprint_npz": str(fingerprint_path.relative_to(ROOT)),
        "events": events,
    }
    (complex_dir / f"{identifier}.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return {
        "payload": payload,
        "bits": bits,
        "features": residue_keys,
        "active_features": binary_features,
        "ifp": ifp,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pH", type=float, default=7.4)
    parser.add_argument("--output-dir", type=Path, default=RESULTS / "prolif" / "corrected")
    parser.add_argument("--pilot-report", type=Path, default=PILOT_REPORT)
    parser.add_argument("--audit-report", type=Path, default=AUDIT_REPORT)
    args = parser.parse_args()
    if not 0.0 <= args.pH <= 14.0:
        parser.error("--pH must be between 0 and 14")
    output_dir = args.output_dir.resolve()
    summary_path = output_dir / "corrected_pilot_summary.json"
    if summary_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite an existing corrected run: {summary_path}"
        )
    for path in (args.pilot_report, args.audit_report, LIGAND_SDF):
        if not path.exists():
            raise FileNotFoundError(path)
    for receptor in RECEPTORS.values():
        for name in ("pdb", "pdbqt", "cif"):
            if not receptor[name].exists():
                raise FileNotFoundError(receptor[name])
    output_dir.mkdir(parents=True, exist_ok=True)
    for subdir in (
        "receptor_preparation", "ligand_reconstruction", "complexes",
        "fingerprints", "matrices", "logs", "reports",
    ):
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)

    baseline = baseline_models(args.pilot_report)
    audited = baseline_models(args.audit_report)
    if set(baseline) != set(audited):
        raise RuntimeError("The audit and original pilot complex selections differ")
    for key in baseline:
        baseline_models_set = {
            model["model"] for model in baseline[key].get("models", [])
        }
        audit_models_set = {
            model["model"] for model in audited[key].get("models", [])
        }
        if baseline_models_set != audit_models_set:
            raise RuntimeError(f"Pose-model selection changed for {key}")
    n_baseline_models = sum(
        len(record.get("models", [])) for record in baseline.values()
    )
    if len(baseline) != 20 or n_baseline_models != 79:
        raise RuntimeError(
            f"Expected the original 20-complex/79-pose pilot, got "
            f"{len(baseline)} complexes/{n_baseline_models} poses"
        )

    sdf_records = read_sdf_records(LIGAND_SDF)
    obabel = shutil.which("obabel")
    if obabel is None:
        raise RuntimeError("Open Babel executable 'obabel' is required for atom mapping")

    receptor_molecules = {}
    receptor_reports = {}
    for receptor_key in RECEPTORS:
        receptor_molecules[receptor_key], receptor_reports[receptor_key] = (
            prepare_receptor(receptor_key, args.pH, output_dir)
        )

    mapping_scratch = output_dir / "logs" / "atom_mapping"
    mapping_scratch.mkdir(parents=True, exist_ok=True)
    unique_ligands = sorted({ligand_id for _, ligand_id in baseline})
    ligand_data = {}
    ligand_reports = {}
    errors = []
    warnings = []
    for ligand_id in unique_ligands:
        if ligand_id not in sdf_records:
            raise FileNotFoundError(f"No curated SDF record for pilot ligand {ligand_id}")
        reference = sdf_records[ligand_id]
        prepared_path = PDBQT_LIGANDS / f"{ligand_id}.pdbqt"
        if not prepared_path.exists():
            raise FileNotFoundError(f"No prepared PDBQT reference for {ligand_id}")
        recovered, mapping, map_report = ligand_atom_mapping(
            reference, prepared_path, obabel, mapping_scratch
        )
        protomer, protomer_mapping, protonation_report = select_protomer(
            reference, args.pH
        )
        heavy_indices = [
            atom.GetIdx() for atom in recovered.GetAtoms()
            if atom.GetSymbol() != "H"
        ]
        prepared_atom_mapping = {
            "ref_to_prepared": mapping,
            "prepared_heavy_indices": heavy_indices,
            "elements": {
                index: recovered.GetAtomWithIdx(index).GetSymbol()
                for index in heavy_indices
            },
        }
        ligand_data[ligand_id] = (
            reference, protomer, protomer_mapping, prepared_atom_mapping
        )
        ligand_reports[ligand_id] = {
            "ligand_id": ligand_id,
            "curated_sdf": str(LIGAND_SDF.relative_to(ROOT)),
            "prepared_docking_pdbqt": str(prepared_path.relative_to(ROOT)),
            "reference_inchikey": reference.GetProp("structure_inchikey")
            if reference.HasProp("structure_inchikey")
            else Chem.MolToInchiKey(Chem.RemoveHs(reference)),
            "original_structure_smiles": reference.GetProp("structure_smiles")
            if reference.HasProp("structure_smiles")
            else Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True),
            "atom_mapping": map_report,
            "protonation": protonation_report,
            "derived_record": f"ligand_reconstruction/{ligand_id}_protonation.json",
        }
        if protonation_report["state_ambiguity"]:
            warnings.append(
                {
                    "ligand_id": ligand_id,
                    "warning": "Dimorphite returned multiple equally ranked protonation states",
                }
            )
        ligand_record_path = (
            output_dir / "ligand_reconstruction" / f"{ligand_id}_protonation.json"
        )
        ligand_record_path.write_text(
            json.dumps(ligand_reports[ligand_id], indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

    category_events = Counter()
    pose_prevalence = Counter()
    ligand_prevalence: dict[str, set[str]] = defaultdict(set)
    all_pose_outputs = []
    matrices: dict[str, list[dict[str, Any]]] = defaultdict(list)
    reproducibility = None
    pose_coordinate_deltas = []
    for (receptor_key, ligand_id), complex_record in sorted(baseline.items()):
        reference, protomer, protomer_mapping, atom_mapping = ligand_data[ligand_id]
        pose_path = Path(complex_record["pose"])
        if not pose_path.is_absolute():
            pose_path = ROOT / pose_path
        if not pose_path.exists():
            errors.append(
                {"receptor_key": receptor_key, "ligand_id": ligand_id, "error": f"missing pose file {pose_path}"}
            )
            continue
        original_atoms = pdbqt_atoms(PDBQT_LIGANDS / f"{ligand_id}.pdbqt")
        pose_models = parse_pose_models(pose_path)
        original_signature = [
            (a["serial"], a["atom_name"], a["element"], a["atom_type"])
            for a in original_atoms
        ]
        for model in complex_record.get("models", []):
            model_number = model["model"]
            matching_model = next(
                (item for item in pose_models if item["model"] == model_number), None
            )
            if matching_model is None:
                errors.append(
                    {
                        "receptor_key": receptor_key,
                        "ligand_id": ligand_id,
                        "model": model_number,
                        "error": "baseline pose model is absent from PDBQT",
                    }
                )
                continue
            pose_signature = [
                (a["serial"], a["atom_name"], a["element"], a["atom_type"])
                for a in matching_model["atoms"]
            ]
            if pose_signature != original_signature:
                errors.append(
                    {
                        "receptor_key": receptor_key,
                        "ligand_id": ligand_id,
                        "model": model_number,
                        "error": "PDBQT atom identities/types/order differ from prepared ligand",
                    }
                )
                continue
            try:
                result = process_one_pose(
                    receptor_key,
                    complex_record,
                    model_number,
                    matching_model["atoms"],
                    reference,
                    protomer,
                    protomer_mapping,
                    atom_mapping,
                    receptor_molecules[receptor_key],
                    output_dir,
                )
                payload = result["payload"]
                payload["pdbqt_input_atom_signature_validated"] = True
                all_pose_outputs.append(payload)
                matrices[receptor_key].append(
                    {
                        "receptor_key": receptor_key,
                        "ligand_id": ligand_id,
                        "pharmacological_label": payload["pharmacological_label"],
                        "pose_model": model_number,
                        "pose_path": str(pose_path.relative_to(ROOT)),
                        "fingerprint_npz": payload["fingerprint_npz"],
                        "bits": result["bits"],
                        "features": result["features"],
                        "active_features": result["active_features"],
                    }
                )
                pose_coordinate_deltas.append(
                    payload["max_heavy_atom_coordinate_delta_angstrom"]
                )
                for interaction, count in payload["interaction_event_counts"].items():
                    category_events[interaction] += count
                    if count:
                        pose_prevalence[interaction] += 1
                        ligand_prevalence[interaction].add(ligand_id)
                if reproducibility is None:
                    check_again = process_one_pose(
                        receptor_key,
                        complex_record,
                        model_number,
                        matching_model["atoms"],
                        reference,
                        protomer,
                        protomer_mapping,
                        atom_mapping,
                        receptor_molecules[receptor_key],
                        output_dir,
                    )
                    reproducibility = {
                        "tested_receptor": receptor_key,
                        "tested_ligand": ligand_id,
                        "tested_pose_model": model_number,
                        "binary_fingerprint_identical": bool(
                            np.array_equal(result["bits"], check_again["bits"])
                        ),
                        "interaction_events_identical": (
                            result["payload"]["events"] == check_again["payload"]["events"]
                        ),
                    }
            except Exception as exc:
                errors.append(
                    {
                        "receptor_key": receptor_key,
                        "ligand_id": ligand_id,
                        "model": model_number,
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                )

    matrix_summaries = {}
    for receptor_key, rows in matrices.items():
        feature_names = sorted(rows[0]["features"]) if rows else []
        matrix = np.zeros((len(rows), len(feature_names)), dtype=np.uint8)
        feature_index = {feature: index for index, feature in enumerate(feature_names)}
        metadata_columns = (
            "receptor_key", "ligand_id", "pharmacological_label", "pose_model",
            "pose_path", "fingerprint_npz",
        )
        csv_path = output_dir / "matrices" / f"IFP_{receptor_key}.csv"
        npz_path = output_dir / "matrices" / f"IFP_{receptor_key}.npz"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow([*metadata_columns, *feature_names])
            for row_index, row in enumerate(rows):
                if sorted(row["features"]) != feature_names:
                    raise RuntimeError(
                        f"Feature schema changed within receptor {receptor_key}"
                    )
                for feature in row["active_features"]:
                    matrix[row_index, feature_index[feature]] = 1
                writer.writerow(
                    [row[column] for column in metadata_columns]
                    + matrix[row_index].tolist()
                )
        np.savez_compressed(
            npz_path,
            matrix=matrix,
            feature_names=np.asarray(feature_names, dtype=str),
            ligand_ids=np.asarray([row["ligand_id"] for row in rows], dtype=str),
            pose_models=np.asarray([row["pose_model"] for row in rows], dtype=np.int16),
        )
        matrix_summaries[receptor_key] = {
            "n_rows": len(rows),
            "n_binary_features": len(feature_names),
            "csv": str(csv_path.relative_to(ROOT)),
            "npz": str(npz_path.relative_to(ROOT)),
        }

    summary_counts = {name: int(category_events[name]) for name in INTERACTIONS}
    baseline_original_counts = baseline_counts(baseline)
    audited_counts = baseline_counts(audited)
    original_by_pose = {
        (receptor, ligand, model["model"]): Counter(
            event.get("interaction", "")
            for event in model.get("interactions", [])
        )
        for (receptor, ligand), record in baseline.items()
        for model in record.get("models", [])
    }
    audit_by_pose = {
        (receptor, ligand, model["model"]): Counter(
            event.get("interaction", "")
            for event in model.get("interactions", [])
        )
        for (receptor, ligand), record in audited.items()
        for model in record.get("models", [])
    }
    docking_input_qc = {}
    for receptor_key in RECEPTORS:
        qc = audit_docking_pdbqt(receptor_key)
        qc["residue_name_mismatch_pose_proximity"] = (
            residue_mismatch_pose_proximity(
                receptor_key,
                qc["residue_name_mismatch_groups"],
                baseline,
            )
        )
        docking_input_qc[receptor_key] = qc
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_complexes_expected": len(baseline),
        "n_complexes_processed": len({(p["receptor_key"], p["ligand_id"]) for p in all_pose_outputs}),
        "n_unique_ligands": len(unique_ligands),
        "n_poses_expected": n_baseline_models,
        "n_poses_processed": len(all_pose_outputs),
        "n_errors": len(errors),
        "n_warnings": len(warnings),
        "interaction_event_counts": summary_counts,
        "poses_with_interaction_category": {
            name: int(pose_prevalence[name]) for name in INTERACTIONS
        },
        "ligands_with_interaction_category": {
            name: len(ligand_prevalence[name]) for name in INTERACTIONS
        },
        "recurrence_definition": (
            "interaction_event_counts count each ProLIF metadata event for a "
            "ligand-residue/protein-residue/type/pose tuple; poses_with_* counts "
            "poses with at least one event of that type; ligands_with_* counts "
            "unique ligand IDs with at least one event within the pilot."
        ),
        "interaction_categories_enabled": list(INTERACTIONS),
        "prolif_version": prolif.__version__,
        "rdkit_version": rdBase.rdkitVersion,
        "dimorphite_version": "2.0.2",
        "pH": args.pH,
        "receptor_preparation": receptor_reports,
        "ligand_preparation": ligand_reports,
        "reproducibility": reproducibility,
        "max_ligand_heavy_coordinate_delta_angstrom": max(pose_coordinate_deltas, default=None),
        "matrix_outputs": matrix_summaries,
        "baseline_comparison": {
            "first_run_values_reported_by_user_not_found_as_independent_artifact": {
                "interaction_event_counts": {
                    "Hydrophobic": 79,
                    "HBDonor": 11,
                    "HBAcceptor": 0,
                    "PiStacking": 0,
                    "Anionic": 0,
                    "Cationic": 0,
                    "CationPi": 0,
                    "PiCation": 0,
                    "VdWContact": 285,
                },
                "total": 375,
                "artifact_available": False,
            },
            "original_pilot": {
                "report": str(args.pilot_report.relative_to(ROOT)),
                "interaction_event_counts": baseline_original_counts,
            },
            "previous_audit": {
                "report": str(args.audit_report.relative_to(ROOT)),
                "interaction_event_counts": audited_counts,
            },
            "corrected_pilot": {"interaction_event_counts": summary_counts},
            "original_vs_audit_same_pose_models": (
                set(original_by_pose) == set(audit_by_pose)
                and all(
                    baseline[key[:2]].get("pose") == audited[key[:2]].get("pose")
                    for key in original_by_pose
                )
            ),
            "original_vs_audit_event_counts_identical": (
                set(original_by_pose) == set(audit_by_pose)
                and all(
                    original_by_pose[key] == audit_by_pose[key]
                    for key in original_by_pose
                )
            ),
        },
        "docking_input_qc": docking_input_qc,
        "ligand_docking_input_qc": {
            ligand_id: audit_ligand_pdbqt(
                PDBQT_LIGANDS / f"{ligand_id}.pdbqt"
            )
            for ligand_id in unique_ligands
        },
        "validation": {
            "same_20_complexes_and_79_model_numbers": (
                len(all_pose_outputs) == 79
                and len({(p["receptor_key"], p["ligand_id"]) for p in all_pose_outputs}) == 20
            ),
            "all_categories_explicitly_enabled": len(INTERACTIONS) == 9,
            "coordinate_preservation": (
                max(pose_coordinate_deltas, default=math.inf) <= 0.001
                and all(
                    receptor["max_heavy_atom_coordinate_delta_angstrom"] <= 0.001
                    for receptor in receptor_reports.values()
                )
            ),
            "reproducibility_identical": bool(
                reproducibility
                and reproducibility["binary_fingerprint_identical"]
                and reproducibility["interaction_events_identical"]
            ),
            "all_poses_processed_without_errors": len(errors) == 0 and len(all_pose_outputs) == 79,
        },
        "decision": "NO_APROBADO",
        "decision_reasons": [
            "La selección de microestado del ligando es una regla determinista sobre variantes Dimorphite y no una población pKa ponderada.",
            "PDB2PQR identificó/añadió átomos pesados ausentes en los PDB de entrada; se descartaron para conservar exactamente el conjunto pesado del receptor usado en docking.",
            "Los PDBQT de 9LL8 y 9LLB reescriben cadenas/numeración y representan algunos residuos parciales con etiquetas distintas; ASP231/ASP232 de 9LLB están a 3.34/4.09 A de una pose piloto, por lo que se recomienda reparar la preparación del receptor y repetir los dockings afectados antes de escalar.",
            "La discrepancia 375/452 entre reportes queda cuantificada, pero no puede atribuirse a una única causa histórica con los artefactos disponibles.",
            "Las frecuencias de interacción no validan por sí solas la química de todos los microestados ni las estructuras incompletas del receptor.",
        ],
        "errors": errors,
        "warnings": warnings,
        "poses": all_pose_outputs,
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    input_paths = [
        LIGAND_SDF,
        args.pilot_report.resolve(),
        args.audit_report.resolve(),
        *[entry[key] for entry in RECEPTORS.values() for key in ("pdb", "pdbqt", "cif")],
        *[
            Path(record["pose"]).resolve()
            for record in baseline.values()
        ],
        Path(__file__).resolve(),
    ]
    inputs = {
        str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path): sha256(path)
        for path in sorted(set(input_paths))
    }
    outputs = {}
    provenance = {
        "generated_at_utc": summary["generated_at_utc"],
        "script": str(Path(__file__).resolve().relative_to(ROOT)),
        "script_sha256": sha256(Path(__file__).resolve()),
        "inputs_sha256": inputs,
        "parameters": {
            "pH": args.pH,
            "interactions": list(INTERACTIONS),
            "prolif_count_fingerprint": False,
            "receptor_forcefield": "AMBER",
            "receptor_pka_method": "PROPKA via PDB2PQR",
            "receptor_heavy_atom_repair_policy": "discard any PDB2PQR-added heavy atom; preserve the input PDB heavy-atom inventory",
            "ligand_protonation_method": "Dimorphite-DL",
            "ligand_state_selection": "max basic-N protonation and carboxyl deprotonation score; deterministic canonical tie-break",
            "coordinate_tolerance_angstrom": 0.001,
        },
        "versions": {
            "python": os.sys.version,
            "rdkit": rdBase.rdkitVersion,
            "prolif": prolif.__version__,
            "dimorphite_dl": "2.0.2",
            "pdb2pqr": subprocess.run(
                [_pdb2pqr_executable(), "--version"],
                capture_output=True,
                text=True,
                check=False,
            ).stdout.strip(),
            "openbabel": subprocess.run(
                [obabel, "-V"], capture_output=True, text=True, check=False
            ).stdout.strip(),
        },
        "derived_outputs_sha256": outputs,
        "manifest_sha256_note": "The manifest excludes its own checksum.",
        "summary_json": str(summary_path.relative_to(ROOT)),
    }
    (output_dir / "provenance_manifest.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    report_lines = [
        "# Corrected ProLIF pilot report",
        "",
        f"**Decision: {summary['decision']} — no procesar todavía los 423 resultados.**",
        "",
        "## Resultado y alcance",
        "",
        f"- Complejos procesados: {summary['n_complexes_processed']}/{summary['n_complexes_expected']}.",
        f"- Modelos de pose procesados: {summary['n_poses_processed']}/{summary['n_poses_expected']}.",
        f"- Errores: {summary['n_errors']}; advertencias registradas: {summary['n_warnings']}.",
        f"- ProLIF: {prolif.__version__}; pH de preparación de receptor/ligando: {args.pH:g}.",
        "- No se repitió docking ni se procesó el conjunto completo.",
        "",
        "## Preparación del receptor",
        "",
        "Se usaron los PDB que alimentaron el docking y PDB2PQR/PROPKA (AMBER, pH indicado, `--nodebump --noopt`). Se conservaron todos los átomos pesados originales y sus coordenadas; átomos pesados reconstruidos por PDB2PQR fueron excluidos para no cambiar el receptor del docking. Se eliminaron enlaces inferidos no peptídicos entre residuos; Asp/Glu se normalizaron según cargas e hidrógenos asignados.",
        "",
        "| Receptor | Átomos pesados fuente | H explícitos | Pesados añadidos por PDB2PQR excluidos | Enlaces inter-residuo espurios retirados | Delta máx. pesada (Å) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for receptor_key, record in receptor_reports.items():
        report_lines.append(
            f"| {receptor_key} | {record['source_heavy_atoms']} | {record['explicit_hydrogens']} | "
            f"{record['pdb2pqr_added_heavy_atoms_discarded']} | "
            f"{len(record['removed_nonlocal_inter_residue_bonds'])} | "
            f"{record['max_heavy_atom_coordinate_delta_angstrom']:.4f} |"
        )
    report_lines.extend(
        [
            "",
            "Residuos diagnósticos (número PDB; cadenas y átomos ausentes quedan explícitos en el resumen JSON):",
            "",
            "| Receptor | Residuo | Presente | Cadena | Átomos pesados faltantes | Átomos aromáticos | Cargas OD1/OD2 ASP155 |",
            "|---|---|---|---|---|---|---|",
        ]
    )
    for receptor_key, record in receptor_reports.items():
        for residue_name, diagnostic in record["diagnostic_residues"].items():
            if diagnostic["present"]:
                missing = ", ".join(diagnostic["missing_expected_heavy_atoms"]) or "ninguno"
                aromatic = ", ".join(diagnostic["aromatic_heavy_atoms"]) or "ninguno"
                charges = diagnostic["carboxyl_oxygen_formal_charges"]
                charge_text = (
                    f"OD1={charges.get('OD1')}, OD2={charges.get('OD2')}"
                    if residue_name == "ASP155"
                    else "—"
                )
                report_lines.append(
                    f"| {receptor_key} | {residue_name} | sí | {diagnostic['chain']} | "
                    f"{missing} | {aromatic} | {charge_text} |"
                )
            else:
                report_lines.append(
                    f"| {receptor_key} | {residue_name} | no | — | — | — | — |"
                )
    report_lines.extend(
        [
            "",
            "La anomalía de 7WC8 se localizó en carboxilatos ASP (incluido ASP217 en la numeración del PDB): RDKit interpretaba una carga O− junto con un doble enlace y reportaba valencia explícita no permitida. Se corrigió el orden de enlace del grupo carboxilato conservando el estado de carga y las coordenadas; los cambios se enumeran por residuo en `receptor_preparation/`. ASP155³·³² se documenta como caso diagnóstico cuando está presente, no como criterio de aceptación.",
            "",
            "PDB2PQR reconstruyó átomos pesados ausentes (que se descartaron) y por ello persiste una limitación: algunos residuos de los PDB fuente pueden estar incompletos. No se minimizó la proteína ni se movieron los átomos pesados.",
            "",
            "## Reconstrucción de ligandos y mapeo",
            "",
            f"Se usaron {len(unique_ligands)} ligandos únicos de `{LIGAND_SDF.relative_to(ROOT)}` para conectividad, órdenes de enlace, aromaticidad y estereoquímica; el PDBQT preparado y cada pose se usaron para identificar átomos y coordenadas. La correspondencia se verificó mediante isomorfismo del grafo pesado y alineación de la conformación de referencia (RMSD registrado por ligando). Los modelos con seriales, tipos o elementos inconsistentes se marcan como error; no se asignan por orden supuesto.",
            "",
            f"Dimorphite-DL enumeró estados a pH {args.pH:g}. Como no entrega poblaciones, la selección aplica una regla documentada de sitios básicos/alifáticos y carboxilos, con desempate canónico; estados alternativos se conservan en el manifiesto por ligando. Esta regla no sustituye un cálculo de población de microestados ni justifica asumir que el docking neutral original es válido para todos los ligandos.",
            "",
            "## Categorías y definición de conteo",
            "",
            "| Categoría | Eventos por pose/residuo/átomo | Poses con ≥1 evento | Ligandos con ≥1 evento |",
            "|---|---:|---:|---:|",
        ]
    )
    for name in INTERACTIONS:
        report_lines.append(
            f"| {name} | {summary_counts[name]} | "
            f"{summary['poses_with_interaction_category'][name]} | "
            f"{summary['ligands_with_interaction_category'][name]} |"
        )
    report_lines.extend(
        [
            "",
            "Un evento es una metadata de interacción devuelta por ProLIF para un par de residuos y una pose; no es una pose única ni un ligando único. La prevalencia de poses y ligandos se informa por separado. Los fingerprints binarios por pose y las matrices IFP se encuentran en `fingerprints/` y `matrices/IFP_<receptor>.*`; cada fila conserva receptor, ligando, etiqueta y número de modelo.",
            "",
            "## Comparación histórica",
            "",
            "| Reporte | Complejos/modelos | Hydrophobic | HBDonor | HBAcceptor | PiStacking | Anionic | Cationic | CationPi | PiCation | VdWContact |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    report_lines.append(
        "| Primera ejecución (valor reportado; artefacto no conservado) | 20/79 | 79 | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 285 |"
    )
    for label, counts in (
        ("JSON de piloto conservado", baseline_original_counts),
        ("Auditoría previa", audited_counts),
        ("Corrected pilot", summary_counts),
    ):
        report_lines.append(
            f"| {label} | 20/79 | "
            + " | ".join(str(counts[name]) for name in INTERACTIONS)
            + " |"
        )
    report_lines.extend(
        [
            "",
            f"En los artefactos JSON actualmente conservados, las rutas/modelos de pose de piloto y auditoría coinciden: `{summary['baseline_comparison']['original_vs_audit_same_pose_models']}`; los eventos por pose también coinciden: `{summary['baseline_comparison']['original_vs_audit_event_counts_identical']}`. Ambos archivos conservados contienen Hydrophobic=156 y total=452. El reporte/archivo de la primera ejecución con Hydrophobic=79 y total=375 no está presente como artefacto independiente con procedencia verificable. Por ello, el aumento de 77 eventos es un dato histórico reportado, pero su causa no se puede atribuir inequívocamente ni descartar mezcla/reescritura de salidas.",
            "",
            "## Reproducibilidad, preparación de Vina y limitaciones",
            "",
            f"Repetición interna del mismo complejo: `{reproducibility}`. Se compara fingerprint binario y eventos serializados. Los hashes de fuentes, poses, reportes y scripts están en `provenance_manifest.json`.",
            "",
            "QC de PDBQT originales de receptor frente a los PDB de docking (correspondencia por coordenadas pesadas; las etiquetas de cadena/número pueden haberse reescrito):",
            "",
            "| Receptor | Coordenadas pesadas coincidentes | Cobertura | Nombres de átomo coincidentes | Nombres de residuo coincidentes | H polares | Cargas no nulas | Renumeración/cadena |",
            "|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for receptor_key, qc in summary["docking_input_qc"].items():
        report_lines.append(
            f"| {receptor_key} | {qc['heavy_atom_coordinate_matches']}/{qc['n_heavy_atoms']} | "
            f"{qc['heavy_atom_coordinate_coverage']:.3f} | "
            f"{qc['atom_name_matches_at_matched_coordinates']}/{qc['heavy_atom_coordinate_matches']} | "
            f"{qc['residue_name_matches_at_matched_coordinates']}/{qc['heavy_atom_coordinate_matches']} | "
            f"{qc['n_polar_hydrogens']} | {qc['nonzero_partial_charge_atoms']} | "
            f"{'sí' if qc['residue_chain_or_numbering_was_rewritten'] else 'no'} |"
        )
        if qc["residue_name_mismatch_groups"]:
            mismatches = ", ".join(
                f"{item['source_residue_name']}{item['source_residue_number']}→"
                f"{item['pdbqt_residue_name']}{item['pdbqt_residue_number']}"
                for item in qc["residue_name_mismatch_groups"]
            )
            closest = qc["residue_name_mismatch_pose_proximity"][:3]
            closest_text = ", ".join(
                f"{item['source_residue_name']}{item['source_residue_number']} "
                f"({item['minimum_distance_to_pilot_pose_angstrom']:.2f} Å; "
                f"{item['nearest_ligand_id']} modelo {item['nearest_pose_model']})"
                for item in closest
            )
            report_lines.append(
                f"Observación {receptor_key}: etiquetas de residuo distintas en las mismas coordenadas: {mismatches}. "
                f"Los tres grupos más próximos a una pose piloto: {closest_text}."
            )
    ligand_qc = list(summary["ligand_docking_input_qc"].values())
    report_lines.extend(
        [
            "",
            f"Los {len(ligand_qc)} PDBQT de ligando inspeccionados tienen en conjunto "
            f"{sum(item['n_polar_hydrogens'] for item in ligand_qc)} hidrógenos polares y "
            f"{sum(item['nonzero_partial_charge_atoms'] for item in ligand_qc)} átomos con carga parcial no nula. "
            "El scoring `vina` estándar no usa cargas parciales electrostáticas, pero sí tipos de átomo y protonación/hidrógenos.",
            "",
            "Los SDF de origen indican `pending_dimorphite_ph7.4`; las aminas neutras pueden cambiar el donador/aceptor y la interacción con ASP155³·³². Los logs revisados identifican Vina 1.2.7 y scoring `vina`; no demuestran que la ejecución fuera Vina-GPU. La conservación de coordenadas no garantiza identidad química completa: los renombrados y residuos incompletos deben verificarse por residuo antes de decidir si redockear.",
            "",
            "El receptor preparado para ProLIF no valida retroactivamente la química de docking. Recomendación: no redockear ahora; primero reconstruir y revisar las cadenas laterales faltantes, especialmente ASP231/ASP232 de 9LLB (3.34/4.09 Å de una pose piloto), y revisar GLU355/LYS223 de 9LL8 (5.61/6.19 Å). Después, repetir los dockings de las combinaciones receptor-ligando cuyas poses queden afectadas. Resolver también los microestados de ligandos básicos: si se confirma que se dockearon neutros y se requiere el estado fisiológico protonado, repetir esos dockings con PDBQT corregidos. No se repitió ningún docking.",
            "",
            "## Archivos de salida y decisión",
            "",
            "- Resumen estructurado: `corrected_pilot_summary.json`.",
            "- Procedencia/hashes: `provenance_manifest.json`.",
            "- Salidas: `receptor_preparation/`, `ligand_reconstruction/`, `complexes/`, `fingerprints/`, `matrices/`, `logs/`.",
            "",
            "**NO APROBADO para procesar los 423 resultados.** La preparación corregida se puede usar como piloto reproducible, pero las alternativas de protonación de ligando, residuos receptor incompletos y ausencia de metadatos exactos de las dos ejecuciones históricas impiden aprobar la interpretación/escalado todavía.",
            "",
        ]
    )
    report_path = output_dir / "reports" / "corrected_pilot_report.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    outputs = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in sorted(output_dir.rglob("*"))
        if path.is_file() and path.name != "provenance_manifest.json"
    }
    provenance["derived_outputs_sha256"] = outputs
    (output_dir / "provenance_manifest.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "n_complexes_processed": summary["n_complexes_processed"],
                "n_poses_processed": summary["n_poses_processed"],
                "n_errors": summary["n_errors"],
                "interaction_event_counts": summary_counts,
                "reproducibility": reproducibility,
                "decision": summary["decision"],
                "report": str(report_path.relative_to(ROOT)),
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
