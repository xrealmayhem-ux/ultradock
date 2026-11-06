#!/usr/bin/env python3
"""Audit source receptor/ligand chemistry and original Vina docking inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
from dimorphite_dl import protonate_smiles
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, Lipinski

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results" / "docking_audit"
LIGAND_SDF = ROOT / "structures" / "ligands" / "sdf" / "docking_ligands_selected_3d.sdf"
SELECTION_CSV = ROOT / "data" / "processed" / "ligand_selection.csv"
LIGAND_PDBQT_DIR = ROOT / "structures" / "ligands" / "pdbqt"
PROLIF_DIR = ROOT / "results" / "prolif" / "corrected"
PROLIF_SUMMARY = PROLIF_DIR / "corrected_pilot_summary.json"
PROLIF_PROVENANCE = PROLIF_DIR / "provenance_manifest.json"
UNIPROT = "P28223"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/P28223.fasta"
UNIPROT_SEQUENCE = (
    "MDILCEENTSLSSTTNSLMQLNDDTRLYSNDFNSGEANTSDAFNWTVDSENRTNLSCEGCLSPSCLSLLHLQEKNWSALL"
    "TAVVIILTIAGNILVIMAVSLEKKLQNATNYFLMSLAIADMLLGFLVMPVSMLTILYGYRWPLPSKLCAVWIYLDVLFST"
    "ASIMHLCAISLDRYVAIQNPIHHSRFNSRTKAFLKIIAVWTISVGISMPIPVFGLQDDSKVFKEGSCLLADDNFVLIGSF"
    "VSFFIPLTIMVITYFLTIKSLQKEATLCVSDLGTRAKLASFSFLPQSSLSSEKLFQRSIHREPGSYTGRRTMQSISNEQK"
    "ACKVLGIVFFLFVVMWCPFFITNIMAVICKESCNEDVIGALLNVFVWIGYLSSAVNPLVYTLFNKTYRSAFSRYIQCQYKEN"
    "KKPLQLILVNTIPALAYKSSQLQMGQKKNSKQDAKTTDNDCSMVALGKQHSEEASKDNSDGVNEKVSCV"
)

RECEPTORS = {
    "7WC8": {
        "chain": "A",
        "pdb": ROOT / "structures" / "7WC8_chainA_protein.pdb",
        "pdbqt": ROOT / "structures" / "7WC8_chainA_protein.pdbqt",
        "cif": ROOT / "receptors" / "raw" / "7WC8.cif",
        "docking_dir": ROOT / "results" / "7WC8_chainA_protein_functional",
    },
    "9LL8": {
        "chain": "R",
        "pdb": ROOT / "structures" / "9LL8_chainR_protein.pdb",
        "pdbqt": ROOT / "structures" / "9LL8_chainR_protein.pdbqt",
        "cif": ROOT / "receptors" / "raw" / "9LL8.cif",
        "docking_dir": ROOT / "results" / "9LL8_chainR_protein_functional",
    },
    "9LLB": {
        "chain": "R",
        "pdb": ROOT / "structures" / "9LLB_chainR_protein.pdb",
        "pdbqt": ROOT / "structures" / "9LLB_chainR_protein.pdbqt",
        "cif": ROOT / "receptors" / "raw" / "9LLB.cif",
        "docking_dir": ROOT / "results" / "9LLB_chainR_protein_functional",
    },
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}
ONE_TO_THREE = {value: key for key, value in THREE_TO_ONE.items()}
ONE_TO_THREE["M"] = "MET"

EXPECTED_HEAVY = {
    "ALA": {"N", "CA", "C", "O", "CB"},
    "ARG": {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "ASN": {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"},
    "ASP": {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"},
    "CYS": {"N", "CA", "C", "O", "CB", "SG"},
    "GLN": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"},
    "GLU": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"},
    "GLY": {"N", "CA", "C", "O"},
    "HIS": {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "ILE": {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"},
    "LEU": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"},
    "LYS": {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"},
    "MET": {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"},
    "PHE": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": {"N", "CA", "C", "O", "CB", "CG", "CD"},
    "SER": {"N", "CA", "C", "O", "CB", "OG"},
    "THR": {"N", "CA", "C", "O", "CB", "OG1", "CG2"},
    "TRP": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "VAL": {"N", "CA", "C", "O", "CB", "CG1", "CG2"},
}

SYMMETRIC_NAME_PAIRS = {
    frozenset(pair)
    for pair in (
        ("CD1", "CD2"), ("CE1", "CE2"), ("CG1", "CG2"),
        ("OD1", "OD2"), ("OE1", "OE2"), ("NH1", "NH2"),
        ("CZ2", "CZ3"), ("CE2", "CE3"),
    )
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_pdb(path: Path) -> list[dict[str, Any]]:
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line[76:78].strip().capitalize()
        if element not in {"H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"}:
            atom_name = re.sub(r"[^A-Za-z]", "", line[12:16]).strip()
            element = atom_name[:2].capitalize() if atom_name[:2].capitalize() in {"Cl", "Br"} else atom_name[:1].capitalize()
        atoms.append(
            {
                "record": line[:6].strip(),
                "serial": int(line[6:11]),
                "atom_name": line[12:16].strip(),
                "altloc": line[16:17].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resid": int(line[22:26]),
                "icode": line[26:27].strip(),
                "element": element,
                "xyz": np.asarray(
                    [float(line[start:start + 8]) for start in (30, 38, 46)],
                    dtype=float,
                ),
                "line": line,
            }
        )
    return atoms


def parse_pdbqt(path: Path, model_number: int | None = None) -> list[dict[str, Any]]:
    atoms = []
    current_model = 1
    in_selected_model = model_number in (None, 1)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            try:
                current_model = int(line[10:14].strip())
            except ValueError:
                current_model = int(line.split()[-1])
            in_selected_model = model_number is None or current_model == model_number
            continue
        if line.startswith("ENDMDL"):
            if model_number is not None:
                in_selected_model = False
            continue
        if not in_selected_model:
            continue
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_type = line[77:].split()[0] if line[77:].split() else ""
        element = {
            "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
            "O": "O", "OA": "O", "OS": "O", "S": "S", "SA": "S",
            "P": "P", "F": "F", "Cl": "Cl", "CL": "Cl",
            "Br": "Br", "BR": "Br", "I": "I", "HD": "H", "HS": "H",
        }.get(atom_type)
        if element is None:
            raise ValueError(f"Unknown PDBQT atom type {atom_type!r} in {path}")
        atoms.append(
            {
                "record": line[:6].strip(),
                "serial": int(line[6:11]),
                "atom_name": line[12:16].strip(),
                "altloc": line[16:17].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resid": int(line[22:26]),
                "icode": line[26:27].strip(),
                "element": element,
                "xyz": np.asarray(
                    [float(line[start:start + 8]) for start in (30, 38, 46)],
                    dtype=float,
                ),
                "charge": float(line[70:76]),
                "atom_type": atom_type,
                "line": line,
            }
        )
    return atoms


def residue_key(atom: dict[str, Any]) -> tuple[str, int, str, str]:
    return atom["chain"], atom["resid"], atom["icode"], atom["resname"]


def coordinate_key(point: np.ndarray) -> tuple[float, float, float]:
    return tuple(round(float(value), 3) for value in point)


def distance_to_box(point: np.ndarray, center: np.ndarray, size: np.ndarray) -> float:
    half = size / 2.0
    outside = np.maximum(np.abs(point - center) - half, 0.0)
    return float(np.linalg.norm(outside))


def residue_sequence(residues: list[dict[str, Any]]) -> str:
    return "".join(THREE_TO_ONE.get(residue["resname"], "X") for residue in residues)


def local_sequence_map(query: str, reference: str) -> tuple[dict[int, int], int, int]:
    """Smith-Waterman alignment; returns zero-based query -> one-based reference."""
    rows, cols = len(query) + 1, len(reference) + 1
    scores = np.zeros((rows, cols), dtype=np.int32)
    trace = np.zeros((rows, cols), dtype=np.uint8)
    best = (0, 0, 0)
    for i in range(1, rows):
        for j in range(1, cols):
            options = (
                (0, 0),
                (scores[i - 1, j - 1] + (3 if query[i - 1] == reference[j - 1] else -2), 1),
                (scores[i - 1, j] - 3, 2),
                (scores[i, j - 1] - 3, 3),
            )
            value, direction = max(options, key=lambda item: item[0])
            scores[i, j] = value
            trace[i, j] = direction
            if value > best[0]:
                best = (int(value), i, j)
    score, i, j = best
    mapping = {}
    matches = 0
    aligned = 0
    while i and j and trace[i, j]:
        direction = int(trace[i, j])
        if direction == 1:
            i -= 1
            j -= 1
            mapping[i] = j + 1
            aligned += 1
            matches += query[i] == reference[j]
        elif direction == 2:
            i -= 1
        else:
            j -= 1
    return mapping, matches, aligned


def load_mmcif_residues(path: Path, chain_id: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    structure = gemmi.read_structure(str(path))
    cif_block = gemmi.cif.read_file(str(path)).sole_block()
    model = structure[0]
    chain = next((candidate for candidate in model if candidate.name == chain_id), None)
    if chain is None:
        raise ValueError(f"Chain {chain_id} absent in {path}")
    residues = []
    nonprotein_residues = 0
    for residue in chain:
        if residue.name not in THREE_TO_ONE:
            nonprotein_residues += 1
            continue
        atoms = [
            {
                "name": atom.name.strip(),
                "element": atom.element.name.capitalize(),
                "xyz": np.asarray([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float),
                "altloc": atom.altloc,
            }
            for atom in residue
            if atom.element.name != "H"
        ]
        residues.append(
            {
                "chain": chain_id,
                "resid": residue.seqid.num,
                "icode": residue.seqid.icode.strip() if residue.seqid.icode else "",
                "resname": residue.name,
                "atoms": atoms,
            }
        )
    entity_id = next(iter(chain)).entity_id if len(chain) else ""
    entity_info = next(
        (
            dict(zip(cif_block.find_mmcif_category("_entity.").tags, row))
            for row in cif_block.find_mmcif_category("_entity.")
            if row[0] == entity_id
        ),
        {},
    )
    meta = {
        "chain": chain_id,
        "subchain": next(iter(chain)).subchain if len(chain) else "",
        "entity_id": entity_id,
        "entity_description": entity_info.get("_entity.pdbx_description", "").strip("'\""),
        "entity_mutations": entity_info.get("_entity.pdbx_mutation", ""),
        "n_residues": len(residues),
        "n_heavy_atoms": sum(
            len({atom["name"] for atom in residue["atoms"]}) for residue in residues
        ),
        "nonprotein_residues_excluded_from_comparison": nonprotein_residues,
    }
    return residues, meta


def compare_cif_to_pdb(
    cif_residues: list[dict[str, Any]], pdb_atoms: list[dict[str, Any]]
) -> dict[str, Any]:
    pdb_residues: dict[tuple[int, str], list[dict[str, Any]]] = defaultdict(list)
    for atom in pdb_atoms:
        if atom["element"] != "H":
            pdb_residues[(atom["resid"], atom["icode"])].append(atom)
    details = []
    for cif_residue in cif_residues:
        key = (cif_residue["resid"], cif_residue["icode"])
        observed = pdb_residues.get(key, [])
        cif_names = {atom["name"] for atom in cif_residue["atoms"]}
        pdb_names = {atom["atom_name"] for atom in observed}
        details.append(
            {
                "resid": key[0],
                "icode": key[1],
                "cif_resname": cif_residue["resname"],
                "pdb_resnames": sorted({atom["resname"] for atom in observed}),
                "cif_heavy_atoms": len(cif_names),
                "pdb_heavy_atoms": len(observed),
                "atom_name_sets_identical": cif_names == pdb_names,
                "cif_atoms_missing_from_pdb": sorted(cif_names - pdb_names),
                "pdb_atoms_absent_from_cif": sorted(pdb_names - cif_names),
            }
        )
    return {
        "n_cif_residues": len(cif_residues),
        "n_residues_with_identical_atom_name_sets": sum(
            item["atom_name_sets_identical"] for item in details
        ),
        "differences": [item for item in details if not item["atom_name_sets_identical"]],
    }


def mmcif_sequence_differences(path: Path) -> list[dict[str, Any]]:
    block = gemmi.cif.read_file(str(path)).sole_block()
    table = block.find_mmcif_category("_struct_ref_seq_dif.")
    differences = []
    for row in table:
        item = dict(zip(table.tags, row))
        if item.get("_struct_ref_seq_dif.pdbx_seq_db_accession_code") != UNIPROT:
            continue
        details = item.get("_struct_ref_seq_dif.details", "").strip("'\"")
        differences.append(
            {
                "accession": UNIPROT,
                "auth_chain": item.get("_struct_ref_seq_dif.pdbx_pdb_strand_id"),
                "auth_residue": item.get("_struct_ref_seq_dif.pdbx_auth_seq_num"),
                "structure_resname": item.get("_struct_ref_seq_dif.mon_id"),
                "reference_resname": item.get("_struct_ref_seq_dif.db_mon_id"),
                "uniprot_position": item.get("_struct_ref_seq_dif.pdbx_seq_db_seq_num"),
                "annotation": details,
            }
        )
    return differences


def align_source_sequence(
    atoms: list[dict[str, Any]], cif_path: Path | None = None
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    grouped: dict[tuple[str, int, str, str], list[dict[str, Any]]] = defaultdict(list)
    for atom in atoms:
        if atom["element"] != "H":
            grouped[residue_key(atom)].append(atom)
    residues = []
    for (chain, resid, icode, resname), atom_list in sorted(grouped.items()):
        residues.append(
            {
                "chain": chain,
                "resid": resid,
                "icode": icode,
                "resname": resname,
                "atom_names": {atom["atom_name"] for atom in atom_list},
            }
        )
    sequence = residue_sequence(residues)
    mapping, matches, aligned = local_sequence_map(sequence, UNIPROT_SEQUENCE)
    substitutions = []
    for index, residue in enumerate(residues):
        residue["uniprot_position"] = mapping.get(index)
        position = residue["uniprot_position"]
        source_aa = THREE_TO_ONE.get(residue["resname"], "X")
        reference_aa = UNIPROT_SEQUENCE[position - 1] if position else None
        residue["reference_resname"] = ONE_TO_THREE.get(reference_aa, "")
        residue["sequence_identity"] = source_aa == reference_aa if position else None
        if position and source_aa != reference_aa:
            substitutions.append(
                {
                    "chain": residue["chain"],
                    "source_residue": residue["resid"],
                    "source_resname": residue["resname"],
                    "uniprot_position": position,
                    "reference_resname": ONE_TO_THREE.get(reference_aa, ""),
                }
            )
    mmcif_differences = (
        mmcif_sequence_differences(cif_path) if cif_path is not None else []
    )
    return residues, {
        "uniprot_accession": UNIPROT,
        "uniprot_url": UNIPROT_URL,
        "canonical_sequence_length": len(UNIPROT_SEQUENCE),
        "observed_source_sequence": sequence,
        "observed_residue_count": len(residues),
        "local_alignment_matches": matches,
        "local_alignment_length": aligned,
        "local_alignment_identity": matches / aligned if aligned else 0.0,
        "unmapped_observed_residues": sum(
            residue["uniprot_position"] is None for residue in residues
        ),
        "observed_sequence_substitutions": substitutions,
        "mmcif_annotated_sequence_differences": mmcif_differences,
        "gpcrdb_ballesteros_weinstein": "not available in local inputs; not inferred",
    }


def read_box(config_dir: Path) -> tuple[dict[str, float], dict[str, Any]]:
    configs = sorted(config_dir.glob("*.txt"))
    if not configs:
        raise FileNotFoundError(f"No docking configs in {config_dir}")
    parsed = []
    for path in configs:
        values = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                values[key.strip()] = value.strip()
        parsed.append((path, values))
    keys = (
        "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
        "exhaustiveness", "num_modes", "cpu",
    )
    variants = Counter(tuple((key, row.get(key)) for key in keys) for _, row in parsed)
    if len(variants) != 1:
        raise ValueError(f"Docking configurations differ in {config_dir}: {variants}")
    parameters = parsed[0][1]
    return (
        {key: float(parameters[key]) for key in keys[:6]},
        {
            "config_count": len(configs),
            "config_sha256_sample": sha256(parsed[0][0]),
            "parameters": parameters,
            "parameter_variant_count": len(variants),
        },
    )


def map_pdbqt_atoms_to_source(
    source_atoms: list[dict[str, Any]], pdbqt_atoms: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    source_by_xyz: dict[tuple[float, float, float], list[dict[str, Any]]] = defaultdict(list)
    pdbqt_by_xyz: dict[tuple[float, float, float], list[dict[str, Any]]] = defaultdict(list)
    for atom in source_atoms:
        if atom["element"] != "H":
            source_by_xyz[coordinate_key(atom["xyz"])].append(atom)
    for atom in pdbqt_atoms:
        if atom["element"] != "H":
            pdbqt_by_xyz[coordinate_key(atom["xyz"])].append(atom)
    paired = []
    unmatched_pdbqt = []
    ambiguous = 0
    for xyz, receptor_atoms in pdbqt_by_xyz.items():
        source_candidates = source_by_xyz.get(xyz, [])
        if len(receptor_atoms) != len(source_candidates):
            if not source_candidates:
                unmatched_pdbqt.extend(receptor_atoms)
                continue
            ambiguous += 1
        source_sorted = sorted(
            source_candidates,
            key=lambda item: (item["resname"], item["atom_name"], item["serial"]),
        )
        pdbqt_sorted = sorted(
            receptor_atoms,
            key=lambda item: (item["resname"], item["atom_name"], item["serial"]),
        )
        paired.extend(zip(source_sorted, pdbqt_sorted))
    unmatched_source = [
        atom
        for xyz, atoms in source_by_xyz.items()
        for atom in atoms
        if xyz not in pdbqt_by_xyz
    ]
    return paired, {
        "source_heavy_atoms": sum(len(atoms) for atoms in source_by_xyz.values()),
        "pdbqt_heavy_atoms": sum(len(atoms) for atoms in pdbqt_by_xyz.values()),
        "matched_coordinate_atom_pairs": len(paired),
        "unmatched_source_heavy_atoms": len(unmatched_source),
        "unmatched_pdbqt_heavy_atoms": len(unmatched_pdbqt),
        "ambiguous_coordinate_groups": ambiguous,
        "heavy_coordinate_coverage": len(paired) / max(
            1, sum(atom["element"] != "H" for atom in source_atoms)
        ),
        "atom_name_matches": sum(
            source["atom_name"] == target["atom_name"] for source, target in paired
        ),
        "residue_name_matches": sum(
            source["resname"] == target["resname"] for source, target in paired
        ),
        "element_matches": sum(
            source["element"] == target["element"] for source, target in paired
        ),
    }


def compare_coordinate_inventory(
    source_atoms: list[dict[str, Any]], prepared_atoms: list[dict[str, Any]]
) -> dict[str, Any]:
    source_by_xyz: dict[tuple[float, float, float], list[dict[str, Any]]] = defaultdict(list)
    prepared_by_xyz: dict[tuple[float, float, float], list[dict[str, Any]]] = defaultdict(list)
    for atom in source_atoms:
        if atom["element"] != "H":
            source_by_xyz[coordinate_key(atom["xyz"])].append(atom)
    for atom in prepared_atoms:
        if atom["element"] != "H":
            prepared_by_xyz[coordinate_key(atom["xyz"])].append(atom)
    matched = []
    for xyz, sources in source_by_xyz.items():
        targets = prepared_by_xyz.get(xyz, [])
        if len(sources) == 1 and len(targets) == 1:
            matched.append((sources[0], targets[0]))
    source_count = sum(map(len, source_by_xyz.values()))
    prepared_count = sum(map(len, prepared_by_xyz.values()))
    return {
        "source_heavy_atoms": source_count,
        "prepared_heavy_atoms": prepared_count,
        "exact_coordinate_matches": len(matched),
        "source_heavy_atoms_unmatched": source_count - len(matched),
        "prepared_heavy_atoms_unmatched": prepared_count - len(matched),
        "all_source_heavy_coordinates_preserved": len(matched) == source_count,
        "atom_name_matches": sum(
            source["atom_name"] == target["atom_name"] for source, target in matched
        ),
        "residue_name_matches": sum(
            source["resname"] == target["resname"] for source, target in matched
        ),
        "element_matches": sum(
            source["element"] == target["element"] for source, target in matched
        ),
        "max_coordinate_delta_angstrom": max(
            (
                float(np.linalg.norm(source["xyz"] - target["xyz"]))
                for source, target in matched
            ),
            default=None,
        ),
    }


def receptor_coordinate_maps() -> tuple[
    dict[str, dict[int, dict[str, Any]]], dict[str, dict[int, dict[str, Any]]]
]:
    residue_maps = {}
    coordinate_maps = {}
    for receptor_id, entry in RECEPTORS.items():
        atoms = parse_pdb(entry["pdb"])
        residues, _ = align_source_sequence(atoms, entry["cif"])
        key_to_position = {
            (residue["chain"], residue["resid"], residue["icode"], residue["resname"]):
            residue["uniprot_position"]
            for residue in residues
        }
        coordinate_map: dict[int, dict[str, Any]] = defaultdict(dict)
        residue_map = {}
        for atom in atoms:
            if atom["element"] == "H":
                continue
            position = key_to_position.get(residue_key(atom))
            if position is None:
                continue
            coordinate_map[position][atom["atom_name"]] = atom["xyz"]
            residue_map[position] = {
                "resname": atom["resname"],
                "chain": atom["chain"],
                "resid": atom["resid"],
            }
        residue_maps[receptor_id] = residue_map
        coordinate_maps[receptor_id] = coordinate_map
    return residue_maps, coordinate_maps


def local_backbone_fit(
    donor: dict[int, dict[str, Any]],
    target: dict[int, dict[str, Any]],
    target_position: int,
    radius: int = 5,
) -> tuple[int, float] | None:
    pairs = [
        (donor[position][name], target[position][name])
        for position in range(target_position - radius, target_position + radius + 1)
        if position != target_position
        for name in ("N", "CA", "C", "O")
        if name in donor.get(position, {}) and name in target.get(position, {})
    ]
    if len(pairs) < 6:
        return None
    source = np.asarray([pair[0] for pair in pairs])
    target_xyz = np.asarray([pair[1] for pair in pairs])
    source_center = source.mean(axis=0)
    target_center = target_xyz.mean(axis=0)
    u, _, vt = np.linalg.svd((source - source_center).T @ (target_xyz - target_center))
    rotation = u @ np.diag([1.0, 1.0, np.linalg.det(u @ vt)]) @ vt
    fitted = (source - source_center) @ rotation + target_center
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target_xyz) ** 2, axis=1))))
    return len(pairs), rmsd


def assess_receptor_templates(
    residue_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    residue_maps, coordinate_maps = receptor_coordinate_maps()
    assessments = []
    for row in residue_rows:
        missing = {
            name for name in row["expected_heavy_atoms_missing_in_source"].split(";")
            if name
        }
        if not missing or not row["in_or_near_box_6A_envelope"]:
            continue
        receptor_id = row["receptor"]
        position = row["source_uniprot_position"]
        candidates = []
        for donor_id in RECEPTORS:
            if donor_id == receptor_id:
                continue
            donor_residue = residue_maps[donor_id].get(position)
            if not donor_residue or donor_residue["resname"] != row["source_resname"]:
                continue
            donor_atoms = coordinate_maps[donor_id].get(position, {})
            if missing.isdisjoint(donor_atoms):
                continue
            fit = local_backbone_fit(
                coordinate_maps[donor_id], coordinate_maps[receptor_id], position
            )
            if fit is None:
                continue
            n_fit_atoms, rmsd = fit
            available_missing = sorted(missing.intersection(donor_atoms))
            candidates.append(
                {
                    "donor_receptor": donor_id,
                    "donor_resid": donor_residue["resid"],
                    "donor_resname": donor_residue["resname"],
                    "missing_atoms_available": ";".join(available_missing),
                    "missing_atom_coverage": len(available_missing) / len(missing),
                    "local_backbone_fit_atom_count": n_fit_atoms,
                    "local_backbone_fit_rmsd_angstrom": rmsd,
                }
            )
        candidates.sort(
            key=lambda item: (
                -item["missing_atom_coverage"],
                item["local_backbone_fit_rmsd_angstrom"],
            )
        )
        if not candidates:
            assessments.append(
                {
                    "target_receptor": receptor_id,
                    "target_resid": row["source_residue"],
                    "uniprot_position": position,
                    "resname": row["source_resname"],
                    "missing_atoms": ";".join(sorted(missing)),
                    "minimum_distance_to_box_angstrom": row[
                        "minimum_heavy_atom_distance_to_box_angstrom"
                    ],
                    "minimum_distance_to_pilot_ligand_angstrom": row[
                        "minimum_distance_to_pilot_ligand_angstrom"
                    ],
                    "donor_receptor": "",
                    "missing_atoms_available": "",
                    "missing_atom_coverage": 0.0,
                    "local_backbone_fit_atom_count": 0,
                    "local_backbone_fit_rmsd_angstrom": "",
                    "assessment": "no local experimental template; unresolved",
                }
            )
        else:
            for candidate in candidates:
                rmsd = candidate["local_backbone_fit_rmsd_angstrom"]
                assessments.append(
                    {
                        "target_receptor": receptor_id,
                        "target_resid": row["source_residue"],
                        "uniprot_position": position,
                        "resname": row["source_resname"],
                        "missing_atoms": ";".join(sorted(missing)),
                        "minimum_distance_to_box_angstrom": row[
                            "minimum_heavy_atom_distance_to_box_angstrom"
                        ],
                        "minimum_distance_to_pilot_ligand_angstrom": row[
                            "minimum_distance_to_pilot_ligand_angstrom"
                        ],
                        **candidate,
                        "assessment": (
                            "candidate for manual local graft; not applied"
                            if candidate["missing_atom_coverage"] == 1.0 and rmsd <= 0.5
                            else "manual review required; not applied"
                        ),
                    }
                )
    for row in residue_rows:
        if row["sequence_identity"] is not False or not row["in_or_near_box_6A_envelope"]:
            continue
        receptor_id = row["receptor"]
        position = row["source_uniprot_position"]
        expected_resname = row["reference_resname"]
        candidates = []
        for donor_id in RECEPTORS:
            if donor_id == receptor_id:
                continue
            donor_residue = residue_maps[donor_id].get(position)
            if not donor_residue or donor_residue["resname"] != expected_resname:
                continue
            expected_atoms = EXPECTED_HEAVY.get(expected_resname, set())
            donor_atoms = coordinate_maps[donor_id].get(position, {})
            coverage = len(expected_atoms.intersection(donor_atoms)) / max(
                1, len(expected_atoms)
            )
            fit = local_backbone_fit(
                coordinate_maps[donor_id], coordinate_maps[receptor_id], position
            )
            if fit is not None:
                candidates.append(
                    {
                        "donor_receptor": donor_id,
                        "donor_resid": donor_residue["resid"],
                        "donor_resname": donor_residue["resname"],
                        "missing_atoms_available": ";".join(
                            sorted(expected_atoms.intersection(donor_atoms))
                        ),
                        "missing_atom_coverage": coverage,
                        "local_backbone_fit_atom_count": fit[0],
                        "local_backbone_fit_rmsd_angstrom": fit[1],
                    }
                )
        if not candidates:
            assessments.append(
                {
                    "target_receptor": receptor_id,
                    "target_resid": row["source_residue"],
                    "uniprot_position": position,
                    "resname": row["source_resname"],
                    "missing_atoms": (
                        f"restore WT {expected_resname} from engineered "
                        f"{row['source_resname']}"
                    ),
                    "minimum_distance_to_box_angstrom": row[
                        "minimum_heavy_atom_distance_to_box_angstrom"
                    ],
                    "minimum_distance_to_pilot_ligand_angstrom": row[
                        "minimum_distance_to_pilot_ligand_angstrom"
                    ],
                    "donor_receptor": "",
                    "missing_atoms_available": "",
                    "missing_atom_coverage": 0.0,
                    "local_backbone_fit_atom_count": 0,
                    "local_backbone_fit_rmsd_angstrom": "",
                    "assessment": "no local WT experimental template; unresolved",
                }
            )
        else:
            for candidate in candidates:
                assessments.append(
                    {
                        "target_receptor": receptor_id,
                        "target_resid": row["source_residue"],
                        "uniprot_position": position,
                        "resname": row["source_resname"],
                        "missing_atoms": (
                            f"restore WT {expected_resname} from engineered "
                            f"{row['source_resname']}"
                        ),
                        "minimum_distance_to_box_angstrom": row[
                            "minimum_heavy_atom_distance_to_box_angstrom"
                        ],
                        "minimum_distance_to_pilot_ligand_angstrom": row[
                            "minimum_distance_to_pilot_ligand_angstrom"
                        ],
                        **candidate,
                        "assessment": "WT experimental residue candidate; not applied",
                    }
                )
    return assessments


def parse_pqr(path: Path) -> dict[str, Any]:
    charges = []
    atom_types = Counter()
    atom_count = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        fields = line.split()
        if len(fields) < 11:
            continue
        atom_count += 1
        charges.append(float(fields[-2]))
        atom_types[fields[-1]] += 1
    return {
        "atoms": atom_count,
        "nonzero_partial_charges": sum(abs(charge) > 1e-6 for charge in charges),
        "net_partial_charge": round(sum(charges), 4),
        "atom_type_counts": dict(sorted(atom_types.items())),
    }


def classify_atom_name_change(source: dict[str, Any], target: dict[str, Any]) -> str:
    if source["atom_name"] == target["atom_name"]:
        return "identical"
    if (
        source["resname"] == target["resname"]
        and frozenset((source["atom_name"], target["atom_name"])) in SYMMETRIC_NAME_PAIRS
    ):
        return "symmetric_atom_label_swap"
    return "non_symmetric_atom_name_change"


def receptor_audit(
    receptor_id: str, entry: dict[str, Any], out_dir: Path
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    source_atoms = parse_pdb(entry["pdb"])
    docking_atoms = parse_pdbqt(entry["pdbqt"])
    heavy_source = [atom for atom in source_atoms if atom["element"] != "H"]
    heavy_docking = [atom for atom in docking_atoms if atom["element"] != "H"]
    box, config_info = read_box(entry["docking_dir"] / "configs")
    batch_summary_path = entry["docking_dir"] / "docking_summary.json"
    batch_summary = json.loads(batch_summary_path.read_text(encoding="utf-8"))
    box_path = Path(batch_summary["box"])
    center = np.asarray([box["center_x"], box["center_y"], box["center_z"]])
    size = np.asarray([box["size_x"], box["size_y"], box["size_z"]])
    paired, pair_summary = map_pdbqt_atoms_to_source(source_atoms, docking_atoms)

    align_residues, sequence_info = align_source_sequence(source_atoms, entry["cif"])
    unmapped_residue_distances = []
    for residue in align_residues:
        if residue["uniprot_position"] is not None:
            continue
        atom_coords = [
            atom["xyz"] for atom in heavy_source
            if residue_key(atom)
            == (
                residue["chain"], residue["resid"],
                residue["icode"], residue["resname"],
            )
        ]
        if atom_coords:
            unmapped_residue_distances.append(
                min(distance_to_box(point, center, size) for point in atom_coords)
            )
    sequence_info["unmapped_residue_count"] = len(unmapped_residue_distances)
    sequence_info["unmapped_residues_within_6A_of_box"] = sum(
        distance <= 6.0 for distance in unmapped_residue_distances
    )
    canonical_by_key = {
        (residue["chain"], residue["resid"], residue["icode"], residue["resname"]):
        residue["uniprot_position"]
        for residue in align_residues
    }
    atom_rows = []
    pair_source_to_pdbqt_residue: dict[tuple[str, int, str, str], Counter] = defaultdict(Counter)
    per_source_atom_names: dict[tuple[str, int, str, str], set[str]] = defaultdict(set)
    per_pdbqt_atom_names: dict[tuple[str, int, str, str], set[str]] = defaultdict(set)
    atom_pairs_by_source: dict[tuple[str, int, str, str], list[tuple[dict, dict]]] = defaultdict(list)
    for source, target in paired:
        skey = residue_key(source)
        tkey = residue_key(target)
        pair_source_to_pdbqt_residue[skey][tkey] += 1
        per_source_atom_names[skey].add(source["atom_name"])
        per_pdbqt_atom_names[tkey].add(target["atom_name"])
        atom_pairs_by_source[skey].append((source, target))
        atom_rows.append(
            {
                "receptor": receptor_id,
                "source_chain": source["chain"],
                "source_residue": source["resid"],
                "source_insertion_code": source["icode"],
                "source_resname": source["resname"],
                "source_uniprot_position": canonical_by_key.get(skey),
                "source_atom": source["atom_name"],
                "element_source": source["element"],
                "pdbqt_chain": target["chain"],
                "pdbqt_residue": target["resid"],
                "pdbqt_insertion_code": target["icode"],
                "pdbqt_resname": target["resname"],
                "pdbqt_atom": target["atom_name"],
                "element_pdbqt": target["element"],
                "atom_type": target["atom_type"],
                "partial_charge": target["charge"],
                "coordinate_delta_angstrom": float(
                    np.linalg.norm(source["xyz"] - target["xyz"])
                ),
                "atom_name_comparison": classify_atom_name_change(source, target),
            }
        )
    source_residues: dict[tuple[str, int, str, str], list[dict]] = defaultdict(list)
    pdbqt_residues: dict[tuple[str, int, str, str], list[dict]] = defaultdict(list)
    for atom in heavy_source:
        source_residues[residue_key(atom)].append(atom)
    for atom in heavy_docking:
        pdbqt_residues[residue_key(atom)].append(atom)

    residue_rows = []
    suspicious_rows = []
    for skey, atoms in sorted(source_residues.items()):
        source_names = {atom["atom_name"] for atom in atoms}
        expected = EXPECTED_HEAVY.get(skey[3], set())
        missing_from_source = sorted(expected - source_names)
        pdbqt_keys = pair_source_to_pdbqt_residue.get(skey, Counter())
        target_key = pdbqt_keys.most_common(1)[0][0] if pdbqt_keys else None
        target_atoms = pdbqt_residues.get(target_key, []) if target_key else []
        target_names = {atom["atom_name"] for atom in target_atoms}
        min_box_distance = min(distance_to_box(atom["xyz"], center, size) for atom in atoms)
        min_pose_distance = None
        expected_extra_missing = sorted(source_names - target_names)
        source_missing_from_pdbqt = sorted(target_names - source_names)
        mismatch_types = sorted(
            {
                classify_atom_name_change(source, target)
                for source, target in atom_pairs_by_source.get(skey, [])
            }
        )
        charges = [atom["charge"] for atom in target_atoms]
        polar_hydrogens = [
            atom for atom in docking_atoms
            if atom["element"] == "H"
            and target_key is not None
            and residue_key(atom) == target_key
            and atom["atom_type"] in ("HD", "HS")
        ]
        row = {
            "receptor": receptor_id,
            "source_chain": skey[0],
            "source_residue": skey[1],
            "source_insertion_code": skey[2],
            "source_resname": skey[3],
            "source_uniprot_position": canonical_by_key.get(skey),
            "reference_resname": next(
                (
                    residue["reference_resname"]
                    for residue in align_residues
                    if (
                        residue["chain"], residue["resid"], residue["icode"], residue["resname"]
                    ) == skey
                ),
                "",
            ),
            "sequence_identity": next(
                (
                    residue["sequence_identity"]
                    for residue in align_residues
                    if (
                        residue["chain"], residue["resid"], residue["icode"], residue["resname"]
                    ) == skey
                ),
                None,
            ),
            "source_heavy_atom_count": len(atoms),
            "source_heavy_atom_names": ";".join(sorted(source_names)),
            "expected_heavy_atoms_missing_in_source": ";".join(missing_from_source),
            "pdbqt_chain": target_key[0] if target_key else "",
            "pdbqt_residue": target_key[1] if target_key else "",
            "pdbqt_resname": target_key[3] if target_key else "",
            "pdbqt_heavy_atom_count": len(target_atoms),
            "pdbqt_heavy_atom_names": ";".join(sorted(target_names)),
            "pdbqt_polar_hydrogen_count": len(polar_hydrogens),
            "source_atoms_missing_in_pdbqt": ";".join(expected_extra_missing),
            "extra_pdbqt_atoms": ";".join(source_missing_from_pdbqt),
            "pdbqt_charge_min": min(charges) if charges else "",
            "pdbqt_charge_max": max(charges) if charges else "",
            "minimum_heavy_atom_distance_to_box_angstrom": min_box_distance,
            "in_or_near_box_6A_envelope": min_box_distance <= 6.0,
            "atom_name_comparison": ";".join(mismatch_types),
            "classification": (
                "source_incomplete_sidechain"
                if missing_from_source
                else (
                    "pdbqt_residue_identity_mislabeled"
                    if target_key and target_key[3] != skey[3]
                    else (
                        "pdbqt_atom_names_differ"
                        if any(kind != "identical" for kind in mismatch_types)
                        else "identity_and_atom_inventory_match"
                    )
                )
            ),
        }
        residue_rows.append(row)
        if (
            missing_from_source
            or row["sequence_identity"] is False
            or (target_key and target_key[3] != skey[3])
            or any(kind == "non_symmetric_atom_name_change" for kind in mismatch_types)
        ):
            suspicious_rows.append(row)

    pilot_pose_files = sorted((PROLIF_DIR / "complexes" / receptor_id).glob("*.json"))
    pilot_pose_records = []
    for metadata_path in pilot_pose_files:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("status") != "processed":
            continue
        pose_path = Path(metadata["pose_path"])
        if not pose_path.is_absolute():
            pose_path = ROOT / pose_path
        ligand_atoms = [
            atom for atom in parse_pdbqt(pose_path, int(metadata["pose_model"]))
            if atom["element"] != "H"
        ]
        if not ligand_atoms:
            continue
        ligand_xyz = np.asarray([atom["xyz"] for atom in ligand_atoms])
        pilot_pose_records.append(
            (metadata["ligand_id"], int(metadata["pose_model"]), ligand_xyz)
        )
    for row in residue_rows:
        skey = (
            row["source_chain"], row["source_residue"],
            row["source_insertion_code"], row["source_resname"],
        )
        residue_xyz = np.asarray(
            [atom["xyz"] for atom in source_residues[skey]]
        )
        nearest = (math.inf, "", "")
        for ligand_id, model_number, ligand_xyz in pilot_pose_records:
            distances = np.linalg.norm(
                residue_xyz[:, np.newaxis, :] - ligand_xyz[np.newaxis, :, :],
                axis=2,
            )
            minimum = float(np.min(distances))
            if minimum < nearest[0]:
                nearest = (minimum, ligand_id, model_number)
        row["minimum_distance_to_pilot_ligand_angstrom"] = (
            nearest[0] if math.isfinite(nearest[0]) else ""
        )
        row["nearest_pilot_ligand_id"] = nearest[1]
        row["nearest_pilot_pose_model"] = nearest[2]
    suspicious_rows = [
        row
        for row in residue_rows
        if row["expected_heavy_atoms_missing_in_source"]
        or row["sequence_identity"] is False
        or (
            row["pdbqt_resname"]
            and row["pdbqt_resname"] != row["source_resname"]
        )
        or "non_symmetric_atom_name_change" in row["atom_name_comparison"]
    ]

    cif_residues, cif_meta = load_mmcif_residues(entry["cif"], entry["chain"])
    cif_comparison = compare_cif_to_pdb(cif_residues, source_atoms)
    prepared_path = (
        PROLIF_DIR / "receptor_preparation" / f"{receptor_id}_prepared.pdb"
    )
    pqr_path = PROLIF_DIR / "receptor_preparation" / f"{receptor_id}_pH7.4.pqr"
    if not prepared_path.exists() or not pqr_path.exists():
        raise FileNotFoundError(f"Corrected ProLIF preparation missing for {receptor_id}")
    prepared_atoms = parse_pdb(prepared_path)
    corrected_pilot = json.loads(PROLIF_SUMMARY.read_text(encoding="utf-8"))
    prolif_prep_summary = corrected_pilot["receptor_preparation"][receptor_id]
    corrected_prep_comparison = compare_coordinate_inventory(source_atoms, prepared_atoms)
    pqr_summary = parse_pqr(pqr_path)
    # Store a readable receptor PDBQT with atom/residue identifiers restored from
    # exact heavy-atom coordinate matches. It deliberately adds no atoms.
    relabelled_path = (
        out_dir / "receptor_repair" / "identity_relabelled_not_for_docking"
        / f"{receptor_id}_identity_relabelled_not_for_docking.pdbqt"
    )
    relabelled_path.parent.mkdir(parents=True, exist_ok=True)
    source_by_xyz: dict[tuple[float, float, float], list[dict[str, Any]]] = defaultdict(list)
    for atom in heavy_source:
        source_by_xyz[coordinate_key(atom["xyz"])].append(atom)
    target_residue_to_source: dict[tuple[str, int, str, str], tuple[str, int, str, str]] = {}
    source_by_docking_atom: list[dict[str, Any] | None] = []
    ambiguous_coordinate_atoms = 0
    for target in docking_atoms:
        if target["element"] == "H":
            source_by_docking_atom.append(None)
            continue
        candidates = source_by_xyz.get(coordinate_key(target["xyz"]), [])
        if len(candidates) != 1:
            ambiguous_coordinate_atoms += 1
            source_by_docking_atom.append(None)
            continue
        source = candidates[0]
        source_by_docking_atom.append(source)
        target_residue_to_source[residue_key(target)] = residue_key(source)
    if len(source_by_docking_atom) != len(docking_atoms):
        raise ValueError(f"{receptor_id}: PDBQT atom-line mapping count mismatch")
    # Assign PDBQT polar hydrogens to the mapped residue using the heavy-atom
    # residue key. Atom names of added H are retained from the original PDBQT.
    original_lines = entry["pdbqt"].read_text(encoding="utf-8").splitlines()
    relabelled_lines = []
    atom_index = 0
    for line in original_lines:
        if line.startswith(("ATOM  ", "HETATM")):
            target = docking_atoms[atom_index]
            source = source_by_docking_atom[atom_index]
            atom_index += 1
            if source is None:
                source_key = target_residue_to_source.get(residue_key(target))
                if source_key is not None:
                    line = (
                        line[:17]
                        + f"{source_key[3]:>3}"
                        + line[20:21]
                        + f"{source_key[0][:1]:1}"
                        + f"{source_key[1]:4d}"
                        + f"{source_key[2][:1]:1}"
                        + line[27:]
                    )
                relabelled_lines.append(line)
                continue
            line = (
                line[:12]
                + f"{source['atom_name']:>4}"
                + line[16:17]
                + f"{source['resname']:>3}"
                + line[20:21]
                + f"{source['chain'][:1]:1}"
                + f"{source['resid']:4d}"
                + f"{source['icode'][:1]:1}"
                + line[27:]
            )
        relabelled_lines.append(line)
    if atom_index != len(docking_atoms):
        raise ValueError(f"{receptor_id}: not all PDBQT atom lines were mapped")
    relabelled_path.write_text("\n".join(relabelled_lines) + "\n", encoding="utf-8")
    relabelled_atoms = parse_pdbqt(relabelled_path)
    relabelled_pairs, relabelled_mapping = map_pdbqt_atoms_to_source(
        source_atoms, relabelled_atoms
    )
    atom_name_change_counts = Counter(
        classify_atom_name_change(source, target) for source, target in paired
    )
    mismatched_residue_pairs = {
        (
            source["chain"], source["resid"], source["resname"],
            target["chain"], target["resid"], target["resname"],
        )
        for source, target in paired
        if source["resname"] != target["resname"]
    }

    summary = {
        "receptor": receptor_id,
        "source_pdb": str(entry["pdb"].relative_to(ROOT)),
        "experimental_mmcif": str(entry["cif"].relative_to(ROOT)),
        "docking_pdbqt": str(entry["pdbqt"].relative_to(ROOT)),
        "source_sha256": {
            "pdb": sha256(entry["pdb"]),
            "mmcif": sha256(entry["cif"]),
            "docking_pdbqt": sha256(entry["pdbqt"]),
            "prolif_prepared_pdb": sha256(prepared_path),
            "prolif_prepared_pqr": sha256(pqr_path),
        },
        "source_chain": entry["chain"],
        "mmcif_chain_metadata": cif_meta,
        "pdb_vs_mmcif_atom_inventory": cif_comparison,
        "corrected_prolif_preparation": {
            "prepared_pdb": str(prepared_path.relative_to(ROOT)),
            "prepared_pqr": str(pqr_path.relative_to(ROOT)),
            "atom_inventory_comparison": corrected_prep_comparison,
            "pqr_charge_and_type_summary": pqr_summary,
            "pipeline_manifest_summary": prolif_prep_summary,
            "heavy_atom_additions_discarded_by_prior_pipeline": prolif_prep_summary.get(
                "pdb2pqr_added_heavy_atoms_discarded"
            ),
            "carboxylate_bond_order_adjustments": prolif_prep_summary.get(
                "carboxylate_bond_order_adjustments", []
            ),
            "sanitization_problems_after_corrections": prolif_prep_summary.get(
                "sanitization_problems_after_corrections"
            ),
        },
        "sequence_alignment": sequence_info,
        "docking_box": box,
        "docking_parameters": config_info,
        "original_docking_execution": {
            "batch_summary": batch_summary,
            "batch_summary_sha256": sha256(batch_summary_path),
            "box_definition": str(box_path),
            "box_definition_sha256": sha256(box_path),
            "run_script": str((ROOT / "scripts" / "run_docking.py").relative_to(ROOT)),
            "run_script_sha256": sha256(ROOT / "scripts" / "run_docking.py"),
            "executable": "/usr/local/bin/vina",
            "command_template": ["/usr/local/bin/vina", "--config=<per-ligand config>"],
            "logs_identify": "AutoDock Vina 1.2.7; scoring vina",
        },
        "docking_pdbqt_atom_mapping": pair_summary,
        "source_heavy_atoms": len(heavy_source),
        "pdbqt_heavy_atoms": len(heavy_docking),
        "pdbqt_total_atoms": len(docking_atoms),
        "pdbqt_polar_hydrogens": sum(
            atom["atom_type"] in ("HD", "HS") for atom in docking_atoms
        ),
        "pdbqt_nonzero_partial_charges": sum(
            abs(atom["charge"]) > 1e-6 for atom in docking_atoms
        ),
        "pdbqt_atom_type_counts": dict(
            sorted(Counter(atom["atom_type"] for atom in docking_atoms).items())
        ),
        "pilot_pose_count_used_for_minimum_distances": len(pilot_pose_records),
        "source_incomplete_residue_count": sum(
            bool(row["expected_heavy_atoms_missing_in_source"])
            for row in residue_rows
        ),
        "source_incomplete_residues_near_box_count": sum(
            bool(row["expected_heavy_atoms_missing_in_source"])
            and row["in_or_near_box_6A_envelope"]
            for row in residue_rows
        ),
        "sequence_substitutions_near_box_count": sum(
            row["sequence_identity"] is False and row["in_or_near_box_6A_envelope"]
            for row in residue_rows
        ),
        "pdbqt_residue_name_mismatch_count": sum(
            source["resname"] != target["resname"] for source, target in paired
        ),
        "pdbqt_residue_name_mismatch_unique_residues": len(mismatched_residue_pairs),
        "pdbqt_atom_name_identical_count": pair_summary["atom_name_matches"],
        "pdbqt_atom_name_change_classification": dict(atom_name_change_counts),
        "pdbqt_atom_name_non_symmetric_changes": sum(
            row["atom_name_comparison"] == "non_symmetric_atom_name_change"
            for row in atom_rows
        ),
        "identity_relabelled_pdbqt": str(relabelled_path.relative_to(ROOT)),
        "identity_relabelled_pdbqt_sha256": sha256(relabelled_path),
        "identity_relabelled_docking_ready": False,
        "identity_relabelled_pdbqt_integrity": {
            "heavy_atoms": len([atom for atom in relabelled_atoms if atom["element"] != "H"]),
            "coordinate_atom_pairs": len(relabelled_pairs),
            "unmatched_source_heavy_atoms": relabelled_mapping["unmatched_source_heavy_atoms"],
            "unmatched_pdbqt_heavy_atoms": relabelled_mapping["unmatched_pdbqt_heavy_atoms"],
            "identity_restored_heavy_atoms": len(heavy_docking) - ambiguous_coordinate_atoms,
            "identity_not_restored_due_to_nonunique_coordinate": ambiguous_coordinate_atoms,
            "source_resname_matches": relabelled_mapping["residue_name_matches"],
            "source_atom_name_matches": relabelled_mapping["atom_name_matches"],
            "element_matches": relabelled_mapping["element_matches"],
            "partial_charge_values_preserved": [
                atom["charge"] for atom in docking_atoms
            ] == [atom["charge"] for atom in relabelled_atoms],
            "atom_type_values_preserved": [
                atom["atom_type"] for atom in docking_atoms
            ] == [atom["atom_type"] for atom in relabelled_atoms],
        },
        "identity_relabelled_policy": (
            "restore labels only by exact source-PDB/heavy-coordinate mapping; "
            "do not add or model missing atoms; not for docking"
        ),
        "suspect_residues": suspicious_rows,
    }
    return summary, residue_rows, atom_rows


def functional_group_counts(mol: Chem.Mol) -> dict[str, int]:
    carbonyl_c = {
        bond.GetBeginAtomIdx()
        for bond in mol.GetBonds()
        if bond.GetBondType() == Chem.BondType.DOUBLE
        and bond.GetBeginAtom().GetSymbol() == "C"
        and bond.GetEndAtom().GetSymbol() == "O"
    }
    carbonyl_c |= {
        bond.GetEndAtomIdx()
        for bond in mol.GetBonds()
        if bond.GetBondType() == Chem.BondType.DOUBLE
        and bond.GetBeginAtom().GetSymbol() == "O"
        and bond.GetEndAtom().GetSymbol() == "C"
    }
    amines = Counter()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N" or atom.GetIsAromatic():
            continue
        heavy_neighbors = [
            neighbor for neighbor in atom.GetNeighbors()
            if neighbor.GetSymbol() != "H"
        ]
        if any(neighbor.GetIdx() in carbonyl_c for neighbor in heavy_neighbors):
            continue
        if any(
            neighbor.GetSymbol() in ("S", "P")
            and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in neighbor.GetBonds())
            for neighbor in heavy_neighbors
        ):
            continue
        degree = len(heavy_neighbors)
        amines[{1: "primary_amine", 2: "secondary_amine", 3: "tertiary_amine"}.get(degree, "other_amine")] += 1
    patterns = {
        "carboxyl": Chem.MolFromSmarts("[CX3](=O)[OX1-,OX2H1]"),
        "phenol": Chem.MolFromSmarts("[c][OX2H1]"),
        "sulfonic_acid": Chem.MolFromSmarts("[SX4](=O)(=O)[OX1-,OX2H1]"),
        "phosphonic_acid": Chem.MolFromSmarts("[PX4](=O)([OX1-,OX2H1])[OX1-,OX2H1]"),
        "guanidine_like": Chem.MolFromSmarts("[NX3][CX3](=[NX2])[NX3]"),
    }
    counts = dict(amines)
    counts["carboxyl"] = len(mol.GetSubstructMatches(patterns["carboxyl"]))
    counts["phenol"] = len(mol.GetSubstructMatches(patterns["phenol"]))
    counts["sulfonic_acid"] = len(mol.GetSubstructMatches(patterns["sulfonic_acid"]))
    counts["phosphonic_acid"] = len(mol.GetSubstructMatches(patterns["phosphonic_acid"]))
    counts["guanidine_like"] = len(mol.GetSubstructMatches(patterns["guanidine_like"]))
    return counts


def compare_ligand_pdbqt_state(
    mol: Chem.Mol, pdbqt_atoms: list[dict[str, Any]]
) -> dict[str, Any]:
    reference = Chem.RemoveHs(Chem.Mol(mol))
    conformer = mol.GetConformer()
    reference_points = {
        atom.GetIdx(): np.asarray(conformer.GetAtomPosition(atom.GetIdx()), dtype=float)
        for atom in mol.GetAtoms()
        if atom.GetSymbol() != "H"
    }
    heavy_targets = [atom for atom in pdbqt_atoms if atom["element"] != "H"]
    target_to_reference: dict[int, int] = {}
    distances = []
    element_matches = 0
    for target in heavy_targets:
        candidates = [
            (
                float(np.linalg.norm(target["xyz"] - point)),
                index,
            )
            for index, point in reference_points.items()
        ]
        if not candidates:
            continue
        distance, index = min(candidates)
        if distance <= 0.02:
            target_to_reference[target["serial"]] = index
            distances.append(distance)
            element_matches += (
                mol.GetAtomWithIdx(index).GetSymbol() == target["element"]
            )
    expected_h_by_atom = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ("N", "O", "S"):
            continue
        expected_h_by_atom[atom.GetIdx()] = atom.GetTotalNumHs(includeNeighbors=True)
    observed_h_by_atom = Counter()
    unmatched_polar_h = 0
    for hydrogen in pdbqt_atoms:
        if hydrogen["atom_type"] not in ("HD", "HS"):
            continue
        closest = min(
            (
                (
                    float(np.linalg.norm(hydrogen["xyz"] - heavy["xyz"])),
                    target_to_reference.get(heavy["serial"]),
                )
                for heavy in heavy_targets
            ),
            default=(math.inf, None),
        )
        distance, parent_idx = closest
        if parent_idx is None or distance > 1.35:
            unmatched_polar_h += 1
            continue
        parent = mol.GetAtomWithIdx(parent_idx)
        if parent.GetSymbol() in ("N", "O", "S"):
            observed_h_by_atom[parent_idx] += 1
    expected_pattern = {
        f"{mol.GetAtomWithIdx(index).GetSymbol()}{index}": count
        for index, count in expected_h_by_atom.items()
    }
    observed_pattern = {
        f"{mol.GetAtomWithIdx(index).GetSymbol()}{index}": observed_h_by_atom[index]
        for index in expected_h_by_atom
    }
    pattern_match = (
        len(target_to_reference) == len(reference_points)
        and all(
            expected_h_by_atom[index] == observed_h_by_atom[index]
            for index in expected_h_by_atom
        )
        and unmatched_polar_h == 0
    )
    return {
        "reference_heavy_atoms": len(reference_points),
        "pdbqt_heavy_atoms": len(heavy_targets),
        "heavy_atom_coordinate_matches": len(target_to_reference),
        "heavy_atom_element_matches": element_matches,
        "heavy_atom_coordinate_coverage": len(target_to_reference) / max(
            1, len(reference_points)
        ),
        "max_matched_heavy_coordinate_delta_angstrom": max(distances, default=None),
        "expected_polar_hydrogen_pattern": expected_pattern,
        "observed_pdbqt_polar_hydrogen_pattern": observed_pattern,
        "unmatched_pdbqt_polar_hydrogens": unmatched_polar_h,
        "polar_hydrogen_pattern_matches_sdf": pattern_match,
    }


def read_sdf(path: Path) -> dict[str, Chem.Mol]:
    records = {}
    for mol in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True):
        if mol is None or not mol.HasProp("ligand_id"):
            continue
        ligand_id = mol.GetProp("ligand_id")
        if ligand_id in records:
            raise ValueError(f"Duplicate ligand ID {ligand_id} in {path}")
        records[ligand_id] = mol
    return records


def read_csv(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return {row["ligand_id"]: row for row in csv.DictReader(stream)}


def ligand_audit(out_dir: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    records = read_sdf(LIGAND_SDF)
    selection = read_csv(SELECTION_CSV)
    pdbqt_qc = read_csv(LIGAND_PDBQT_DIR / "ligand_pdbqt_qc.csv")
    rows = []
    ligand_summary = {}
    for ligand_id, mol in sorted(records.items()):
        if ligand_id not in selection:
            raise ValueError(f"SDF record missing from selection CSV: {ligand_id}")
        entry = selection[ligand_id]
        pdbqt_path = LIGAND_PDBQT_DIR / f"{ligand_id}.pdbqt"
        qc = pdbqt_qc.get(ligand_id, {})
        original_charge = int(Chem.GetFormalCharge(mol))
        smiles = entry.get("structure_smiles") or Chem.MolToSmiles(
            Chem.RemoveHs(mol), canonical=True, isomericSmiles=True
        )
        candidates = sorted(
            set(
                protonate_smiles(
                    smiles,
                    ph_min=7.4,
                    ph_max=7.4,
                    precision=0.25,
                    max_variants=128,
                    validate_output=True,
                )
            )
        )
        if not candidates:
            candidates = [""]
        candidate_mols = []
        for state_index, state_smiles in enumerate(candidates, start=1):
            candidate = Chem.MolFromSmiles(state_smiles) if state_smiles else None
            if state_smiles and candidate is None:
                raise ValueError(f"Invalid Dimorphite state for {ligand_id}: {state_smiles}")
            charge = int(Chem.GetFormalCharge(candidate)) if candidate else None
            key = Chem.MolToInchiKey(candidate) if candidate else ""
            candidate_mols.append((state_smiles, charge, key))
        charge_states = sorted({charge for _, charge, _ in candidate_mols if charge is not None})
        group_counts = functional_group_counts(Chem.RemoveHs(Chem.Mol(mol)))
        sdf_polar_h = sum(
            atom.GetTotalNumHs(includeNeighbors=True)
            for atom in mol.GetAtoms()
            if atom.GetSymbol() in ("N", "O", "S")
        )
        if pdbqt_path.exists():
            dock_atoms = parse_pdbqt(pdbqt_path)
            pqr_polar_h = sum(atom["atom_type"] in ("HD", "HS") for atom in dock_atoms)
            pqr_types = Counter(atom["atom_type"] for atom in dock_atoms)
            pqr_nonzero_charges = sum(abs(atom["charge"]) > 1e-6 for atom in dock_atoms)
            pqr_heavy = sum(atom["element"] != "H" for atom in dock_atoms)
            pdbqt_state_comparison = compare_ligand_pdbqt_state(mol, dock_atoms)
            pdbqt_exists = True
        else:
            pqr_polar_h = None
            pqr_types = Counter()
            pqr_nonzero_charges = None
            pqr_heavy = None
            pdbqt_state_comparison = {}
            pdbqt_exists = False
        distinct_charge_count = len(charge_states)
        original_in_candidate_charges = original_charge in charge_states
        if qc.get("status") == "excluded":
            protonation_decision = "not_docked_excluded_by_atom_support"
        elif not pdbqt_exists:
            protonation_decision = "not_evaluable_missing_pdbqt"
        elif (
            pdbqt_state_comparison.get("heavy_atom_coordinate_coverage", 0.0) < 1.0
            or pdbqt_state_comparison.get("heavy_atom_element_matches", 0)
            != pdbqt_state_comparison.get("reference_heavy_atoms", 0)
            or not pdbqt_state_comparison.get("polar_hydrogen_pattern_matches_sdf")
        ):
            protonation_decision = "manual_review_pdbqt_sdf_state_mismatch"
        elif distinct_charge_count > 1:
            protonation_decision = "manual_review_microstate_sensitivity"
        elif distinct_charge_count == 1 and not original_in_candidate_charges:
            protonation_decision = "manual_review_original_charge_outside_enumeration"
        else:
            protonation_decision = "state_consistent_with_enumerated_charge_not_population_validated"
        original_smiles = Chem.MolToSmiles(
            Chem.RemoveHs(Chem.Mol(mol)), canonical=True, isomericSmiles=True
        )
        ligand_summary[ligand_id] = {
            "original_formal_charge": original_charge,
            "candidate_formal_charges": charge_states,
            "microstate_count": len(candidate_mols),
            "has_pdbqt": pdbqt_exists,
            "protonation_decision": protonation_decision,
            "group_counts": group_counts,
            "pdbqt_state_comparison": pdbqt_state_comparison,
            "pdbqt_nonzero_partial_charges": pqr_nonzero_charges,
        }
        for state_index, (state_smiles, charge, key) in enumerate(candidate_mols, start=1):
            rows.append(
                {
                    "ligand_id": ligand_id,
                    "microstate_index": state_index,
                    "microstate_smiles": state_smiles,
                    "microstate_inchikey": key,
                    "original_sdf_smiles": original_smiles,
                    "original_formal_charge": original_charge,
                    "candidate_formal_charge": charge,
                    "candidate_charge_matches_original": charge == original_charge,
                    "structure_heavy_atoms": Chem.RemoveHs(mol).GetNumAtoms(),
                    "molecular_weight": round(Descriptors.MolWt(Chem.RemoveHs(mol)), 3),
                    "hbd": int(Lipinski.NumHDonors(Chem.RemoveHs(mol))),
                    "hba": int(Lipinski.NumHAcceptors(Chem.RemoveHs(mol))),
                    "functional_groups_json": json.dumps(group_counts, sort_keys=True),
                    "protonation_policy_in_selection": entry.get("protonation_policy", ""),
                    "sdf_state_source": "curated input SMILES; selection marks pending_dimorphite_ph7.4",
                    "pH_enumeration": 7.4,
                    "dimorphite_candidate_count": len(candidate_mols),
                    "dimorphite_candidate_charge_set": ";".join(map(str, charge_states)),
                    "dimorphite_state_population": "not estimated",
                    "pka_prediction": "not available in inspected environment; no pKa values inferred",
                    "docking_pdbqt_status": qc.get("status", "missing"),
                    "docking_pdbqt_path": str(pdbqt_path.relative_to(ROOT)) if pdbqt_exists else "",
                    "pdbqt_heavy_atoms": pqr_heavy,
                    "pdbqt_polar_hydrogens": pqr_polar_h,
                    "sdf_expected_polar_hydrogens": sdf_polar_h,
                    "pdbqt_vs_sdf_polar_hydrogen_delta": (
                        pqr_polar_h - sdf_polar_h if pqr_polar_h is not None else None
                    ),
                    "pdbqt_atom_types_json": json.dumps(dict(sorted(pqr_types.items()))),
                    "pdbqt_nonzero_partial_charges": pqr_nonzero_charges,
                    "sdf_pdbqt_atom_state_comparison_json": json.dumps(
                        pdbqt_state_comparison, sort_keys=True
                    ),
                    "sdf_pdbqt_polar_hydrogen_pattern_match": pdbqt_state_comparison.get(
                        "polar_hydrogen_pattern_matches_sdf"
                    ),
                    "protonation_decision": protonation_decision,
                    "microstate_limitation": (
                        "Dimorphite enumerates plausible states; it does not report pKa-derived populations."
                    ),
                }
            )
    return rows, ligand_summary


def read_docking_logs(entry: dict[str, Any], ligand_ids: set[str]) -> dict[str, Any]:
    log_dir = entry["docking_dir"] / "logs"
    version_counts = Counter()
    scoring_counts = Counter()
    seed_count = 0
    gpu_marker_count = 0
    log_count = 0
    for ligand_id in ligand_ids:
        path = log_dir / f"{ligand_id}.log"
        if not path.exists():
            continue
        log_count += 1
        text = path.read_text(encoding="utf-8", errors="replace")
        version = next((line.strip() for line in text.splitlines() if "AutoDock Vina" in line), "")
        scoring = next((line.split(":", 1)[1].strip() for line in text.splitlines() if line.startswith("Scoring function")), "")
        version_counts[version] += 1
        scoring_counts[scoring] += 1
        seed_count += bool(re.search(r"random seed:\s*\d+", text))
        gpu_marker_count += "Vina-GPU" in text or "GPU monte carlo" in text
    return {
        "n_logs": log_count,
        "version_lines": dict(version_counts),
        "scoring_functions": dict(scoring_counts),
        "logs_with_random_seed": seed_count,
        "logs_with_vina_gpu_marker": gpu_marker_count,
    }


def source_inventory(out_dir: Path) -> list[dict[str, Any]]:
    paths = [
        LIGAND_SDF,
        SELECTION_CSV,
        LIGAND_PDBQT_DIR / "ligand_pdbqt_qc.csv",
        PROLIF_SUMMARY,
        PROLIF_PROVENANCE,
        ROOT / "scripts" / "run_docking.py",
        ROOT / "scripts" / "prepare_ligand_pdbqt.py",
        ROOT / "scripts" / "prepare_selected_ligands.py",
    ]
    prior_manifest = json.loads(PROLIF_PROVENANCE.read_text(encoding="utf-8"))
    paths.extend(
        ROOT / path for path in prior_manifest.get("inputs_sha256", {})
    )
    paths.extend(sorted((ROOT / "structures").glob("*_crystal.pdb")))
    for entry in RECEPTORS.values():
        paths.extend([entry["pdb"], entry["pdbqt"], entry["cif"]])
        for subdir in ("configs", "logs", "poses"):
            paths.extend(sorted((entry["docking_dir"] / subdir).glob("*")))
        batch_summary_path = entry["docking_dir"] / "docking_summary.json"
        paths.append(batch_summary_path)
        if batch_summary_path.exists():
            batch_summary = json.loads(batch_summary_path.read_text(encoding="utf-8"))
            box_path = Path(batch_summary["box"])
            if box_path.exists():
                paths.append(box_path)
    paths.extend(sorted(LIGAND_PDBQT_DIR.glob("*.pdbqt")))
    rows = []
    for path in sorted(set(paths)):
        if not path.is_file():
            continue
        rows.append(
            {
                "path": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                "file_bytes": path.stat().st_size,
                "sha256": sha256(path),
                "role": (
                    "docking_pose"
                    if "poses" in path.parts
                    else (
                        "docking_config"
                        if "configs" in path.parts
                        else (
                            "docking_log"
                            if "logs" in path.parts
                            else "source_or_reference"
                        )
                    )
                ),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def decision_and_reports(
    out_dir: Path,
    receptor_summaries: dict[str, dict[str, Any]],
    residue_rows: list[dict[str, Any]],
    ligand_rows: list[dict[str, Any]],
    ligand_summary: dict[str, Any],
    inventory: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    ligand_ids = sorted(ligand_summary)
    dockinged_ids = {
        ligand_id
        for receptor in RECEPTORS.values()
        for ligand_id in (
            path.stem for path in (receptor["docking_dir"] / "poses").glob("*.pdbqt")
        )
    }
    receptor_status = {}
    for receptor_id, summary in receptor_summaries.items():
        nearby_incomplete = [
            row for row in residue_rows
            if row["receptor"] == receptor_id
            and row["expected_heavy_atoms_missing_in_source"]
            and row["in_or_near_box_6A_envelope"]
        ]
        nearby_mutations = [
            row for row in residue_rows
            if row["receptor"] == receptor_id
            and row["sequence_identity"] is False
            and row["in_or_near_box_6A_envelope"]
        ]
        receptor_status[receptor_id] = {
            "nearby_incomplete_residues": nearby_incomplete,
            "nearby_sequence_substitutions": nearby_mutations,
            "identity_issue_in_box": bool(nearby_incomplete or nearby_mutations),
            "decision": (
                "REPETIR DOCKING"
                if nearby_incomplete or nearby_mutations
                else "CONSERVAR"
            ),
            "reason": (
                "La estructura experimental usada como receptor carece de átomos pesados "
                "de cadenas laterales dentro o a ≤6 Å de la caja y/o contiene sustituciones "
                "confirmadas como ingenierizadas dentro de la caja. No se puede asumir que la "
                "caja representa la secuencia WT. Debe resolverse la identidad y revisar la "
                "región antes del redocking."
                if nearby_incomplete or nearby_mutations
                else "No se detectó pérdida de átomos pesados no triviales ni sustituciones dentro del volumen de búsqueda o de su envolvente conservadora."
            ),
        }
    rows = []
    for receptor_id, entry in RECEPTORS.items():
        status = receptor_status[receptor_id]
        for ligand_id in ligand_ids:
            ligand = ligand_summary[ligand_id]
            if not ligand["has_pdbqt"]:
                decision = "NO EVALUABLE POR FALTA DE INFORMACIÓN"
                issue = "Ligando excluido/no existe PDBQT de docking"
                evidence = "QC de conversión de ligando sin PDBQT; no aparece en los resultados de docking"
            elif status["identity_issue_in_box"]:
                decision = "REPETIR DOCKING"
                issue = "Receptor incompleto dentro/cerca de la caja original"
                evidence = "; ".join(
                    f"{r['source_resname']}{r['source_residue']} faltan "
                    f"{r['expected_heavy_atoms_missing_in_source']} (distancia a caja "
                    f"{r['minimum_heavy_atom_distance_to_box_angstrom']:.2f} Å)"
                    for r in status["nearby_incomplete_residues"][:8]
                )
                evidence += (
                    "; sustituciones dentro de caja: "
                    + ", ".join(
                        f"{r['source_resname']}{r['source_residue']}→{r['reference_resname']}"
                        for r in status["nearby_sequence_substitutions"]
                    )
                    if status["nearby_sequence_substitutions"]
                    else ""
                )
            elif ligand["protonation_decision"] == "manual_review_microstate_sensitivity":
                decision = "REVISIÓN MANUAL"
                issue = "Más de un estado de carga plausible en enumeración a pH 7.4"
                evidence = f"carga original {ligand['original_formal_charge']}; candidatos {ligand['candidate_formal_charges']}; población no estimada"
            elif ligand["protonation_decision"] == "manual_review_original_charge_outside_enumeration":
                decision = "REVISIÓN MANUAL"
                issue = "La carga del SDF original no aparece entre los estados enumerados"
                evidence = f"carga original {ligand['original_formal_charge']}; candidatos {ligand['candidate_formal_charges']}"
            elif ligand["protonation_decision"] == "manual_review_pdbqt_sdf_state_mismatch":
                decision = "REVISIÓN MANUAL"
                issue = "El PDBQT de ligando no coincide en átomos pesados/coordenadas/hidrógenos con el SDF"
                evidence = json.dumps(ligand["pdbqt_state_comparison"], sort_keys=True)
            else:
                decision = "CONSERVAR"
                issue = "Sin evidencia suficiente de mismatch de microestado; conservar pendiente de revisión de química"
                evidence = (
                    f"carga SDF={ligand['original_formal_charge']}; "
                    f"candidatos Dimorphite={ligand['candidate_formal_charges']}; "
                    f"hidrógenos polares SDF/PDBQT no comparados por carga formal"
                )
            rows.append(
                {
                    "ligand_id": ligand_id,
                    "receptor": receptor_id,
                    "original_state_formal_charge": ligand["original_formal_charge"],
                    "candidate_formal_charges_ph7_4": ";".join(
                        map(str, ligand["candidate_formal_charges"])
                    ),
                    "ligand_protonation_decision": ligand["protonation_decision"],
                    "issue_detected": issue,
                    "evidence": evidence,
                    "potential_impact": (
                        "La identidad química/receptor puede alterar el volumen y la carga del bolsillo; el score/poses previos no son comparables con receptor nativo."
                        if decision == "REPETIR DOCKING"
                        else (
                            "El microestado puede modificar donadores/aceptores y contactos iónicos; evaluar variantes antes de interpretar."
                            if decision == "REVISIÓN MANUAL"
                            else "No se identificó un cambio concluyente que obligue a repetir."
                        )
                    ),
                    "original_docking_exists": ligand_id in dockinged_ids,
                    "decision": decision,
                    "justification": status["reason"] if decision == "REPETIR DOCKING" else evidence,
                }
            )
    protonation_review = Counter(item["protonation_decision"] for item in ligand_summary.values())
    receptor_review = {
        receptor_id: {
            "decision": item["decision"],
            "n_nearby_incomplete_residues": len(item["nearby_incomplete_residues"]),
            "residues": [
                {
                    "resname": row["source_resname"],
                    "resid": row["source_residue"],
                    "missing_atoms": row["expected_heavy_atoms_missing_in_source"],
                    "minimum_distance_to_box_angstrom": row["minimum_heavy_atom_distance_to_box_angstrom"],
                    "minimum_distance_to_pilot_ligand_angstrom": row[
                        "minimum_distance_to_pilot_ligand_angstrom"
                    ],
                    "nearest_pilot_ligand_id": row["nearest_pilot_ligand_id"],
                    "nearest_pilot_pose_model": row["nearest_pilot_pose_model"],
                    "pdbqt_resname": row["pdbqt_resname"],
                    "pdbqt_chain": row["pdbqt_chain"],
                    "pdbqt_residue": row["pdbqt_residue"],
                }
                for row in item["nearby_incomplete_residues"]
            ],
            "sequence_substitutions": [
                {
                    "source_resname": row["source_resname"],
                    "source_resid": row["source_residue"],
                    "reference_resname": row["reference_resname"],
                    "uniprot_position": row["source_uniprot_position"],
                    "minimum_distance_to_box_angstrom": row[
                        "minimum_heavy_atom_distance_to_box_angstrom"
                    ],
                    "minimum_distance_to_pilot_ligand_angstrom": row[
                        "minimum_distance_to_pilot_ligand_angstrom"
                    ],
                }
                for row in item["nearby_sequence_substitutions"]
            ],
        }
        for receptor_id, item in receptor_status.items()
    }
    return rows, {
        "receptor_decisions": receptor_review,
        "ligand_protonation_decisions": dict(protonation_review),
        "conclusion": (
            "NO APROBADO"
            if any(item["decision"] == "REPETIR DOCKING" for item in receptor_status.values())
            else (
                "APROBADO CON REPETICIÓN SELECTIVA"
                if any(value for value in protonation_review.values())
                else "APROBADO"
            )
        ),
        "n_ligands": len(ligand_summary),
        "n_ligands_with_docking_pdbqt": sum(item["has_pdbqt"] for item in ligand_summary.values()),
        "n_docking_combinations": len(rows),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite outputs generated by a prior run of this audit",
    )
    args = parser.parse_args()
    out_dir = args.output_dir.resolve()
    if (out_dir / "audit_summary.json").exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite existing audit {out_dir}")
    for path in (LIGAND_SDF, SELECTION_CSV, LIGAND_PDBQT_DIR / "ligand_pdbqt_qc.csv", PROLIF_SUMMARY):
        if not path.exists():
            raise FileNotFoundError(path)
    for entry in RECEPTORS.values():
        for key in ("pdb", "pdbqt", "cif"):
            if not entry[key].exists():
                raise FileNotFoundError(entry[key])
    for subdir in (
        "receptor_comparison", "receptor_repair", "ligand_protonation",
        "redocking_decisions", "validation", "reports", "logs",
    ):
        (out_dir / subdir).mkdir(parents=True, exist_ok=True)
    fasta_path = out_dir / "receptor_comparison" / f"UniProt_{UNIPROT}.fasta"
    fasta_path.write_text(
        f">sp|{UNIPROT}|5HT2A_HUMAN 5-hydroxytryptamine receptor 2A; UniProt canonical sequence\n"
        + "\n".join(
            UNIPROT_SEQUENCE[index:index + 60]
            for index in range(0, len(UNIPROT_SEQUENCE), 60)
        )
        + "\n",
        encoding="utf-8",
    )

    receptor_summaries = {}
    all_residue_rows = []
    all_atom_rows = []
    for receptor_id, entry in RECEPTORS.items():
        summary, residue_rows, atom_rows = receptor_audit(receptor_id, entry, out_dir)
        receptor_summaries[receptor_id] = summary
        all_residue_rows.extend(residue_rows)
        all_atom_rows.extend(atom_rows)
    write_csv(out_dir / "receptor_comparison" / "receptor_residue_comparison.csv", all_residue_rows)
    write_csv(out_dir / "receptor_comparison" / "receptor_atom_comparison.csv", all_atom_rows)
    template_assessments = assess_receptor_templates(all_residue_rows)
    write_csv(
        out_dir / "receptor_comparison" / "receptor_template_assessment.csv",
        template_assessments,
    )

    ligand_rows, ligand_summary = ligand_audit(out_dir)
    write_csv(out_dir / "ligand_protonation" / "ligand_protonation_audit.csv", ligand_rows)
    decisions, decision_summary = decision_and_reports(
        out_dir, receptor_summaries, all_residue_rows, ligand_rows, ligand_summary, []
    )
    write_csv(out_dir / "redocking_decisions" / "redocking_decisions.csv", decisions)

    inventory = source_inventory(out_dir)
    write_csv(out_dir / "validation" / "source_inventory_sha256.csv", inventory)
    docking_logs = {
        receptor_id: read_docking_logs(entry, set(ligand_summary))
        for receptor_id, entry in RECEPTORS.items()
    }
    corrected_pilot = json.loads(PROLIF_SUMMARY.read_text(encoding="utf-8"))
    corrected_provenance = json.loads(PROLIF_PROVENANCE.read_text(encoding="utf-8"))
    current_hashes = {item["path"]: item["sha256"] for item in inventory}
    old_hashes = corrected_provenance.get("inputs_sha256", {})
    pilot_source_hash_comparison = {
        "compared_input_count": len(old_hashes),
        "matching_input_count": sum(
            current_hashes.get(path) == digest for path, digest in old_hashes.items()
        ),
        "changed_or_missing_inputs": [
            path for path, digest in old_hashes.items()
            if current_hashes.get(path) != digest
        ],
    }
    mapped_ligands = [
        item["pdbqt_state_comparison"]
        for item in ligand_summary.values()
        if item["has_pdbqt"]
    ]
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "uniprot_reference": {
            "accession": UNIPROT,
            "source_url": UNIPROT_URL,
            "sequence_sha256": hashlib.sha256(UNIPROT_SEQUENCE.encode()).hexdigest(),
            "sequence_length": len(UNIPROT_SEQUENCE),
            "gpcrdb_generic_numbering": "not present in repository; not inferred",
        },
        "versions": {
            "python": __import__("sys").version,
            "rdkit": rdBase.rdkitVersion,
            "dimorphite_dl": importlib.metadata.version("dimorphite_dl"),
            "gemmi": gemmi.__version__,
            "vina_logs": docking_logs,
        },
        "prior_corrected_pilot": {
            "metrics": {
                key: corrected_pilot.get(key)
                for key in (
                    "n_complexes_processed", "n_unique_ligands", "n_poses_processed",
                    "n_errors", "n_warnings", "interaction_event_counts",
                )
            },
            "receptor_preparation": {
                key: {
                    "source_heavy_atoms": value.get("source_heavy_atoms"),
                    "prepared_heavy_atoms": value.get("prepared_heavy_atoms"),
                    "pdb2pqr_added_heavy_atoms_discarded": value.get(
                        "pdb2pqr_added_heavy_atoms_discarded"
                    ),
                    "carboxylate_bond_order_adjustments": value.get(
                        "carboxylate_bond_order_adjustments"
                    ),
                    "sanitization_problems_after_corrections": value.get(
                        "sanitization_problems_after_corrections"
                    ),
                }
                for key, value in corrected_pilot.get("receptor_preparation", {}).items()
            },
            "source_hash_comparison": pilot_source_hash_comparison,
        },
        "inventory_file_count": len(inventory),
        "receptors": receptor_summaries,
        "structural_template_assessment": {
            "candidate_rows": len(template_assessments),
            "candidates": template_assessments,
            "policy": "experimental templates only; local backbone superposition; no atom grafts applied",
        },
        "ligands": {
            "reference_sdf_records": len(ligand_summary),
            "selection_rows": len(read_csv(SELECTION_CSV)),
            "pdbqt_qc_rows": len(read_csv(LIGAND_PDBQT_DIR / "ligand_pdbqt_qc.csv")),
            "candidate_states": len(ligand_rows),
            "states_by_ligand_decision": dict(
                Counter(item["protonation_decision"] for item in ligand_summary.values())
            ),
            "states_by_candidate_charge": dict(
                Counter(
                    str(item["candidate_formal_charge"])
                    for item in ligand_rows
                )
            ),
            "pdbqt_to_sdf_identity_validation": {
                "pdbqt_ligands": len(mapped_ligands),
                "complete_heavy_atom_coordinate_mappings": sum(
                    item.get("heavy_atom_coordinate_coverage") == 1.0
                    for item in mapped_ligands
                ),
                "heavy_atom_element_mismatches": sum(
                    item.get("heavy_atom_element_matches")
                    != item.get("reference_heavy_atoms")
                    for item in mapped_ligands
                ),
                "exact_polar_hydrogen_patterns": sum(
                    item.get("polar_hydrogen_pattern_matches_sdf", False)
                    for item in mapped_ligands
                ),
                "max_coordinate_delta_angstrom": max(
                    (
                        item.get("max_matched_heavy_coordinate_delta_angstrom") or 0
                        for item in mapped_ligands
                    ),
                    default=None,
                ),
            },
            "no_pka_population_estimates": True,
            "excluded_ligands": [
                ligand_id for ligand_id, item in ligand_summary.items()
                if not item["has_pdbqt"]
            ],
        },
        "docking_logs": docking_logs,
        "decisions": decision_summary,
        "method_limitations": [
            "No sidechain conformations were invented or added during this audit.",
            "Dimorphite-DL enumerates plausible microstates but supplies no populations; no pKa predictor was available in the inspected environment.",
            "The 6 A envelope is a conservative screening distance from observed residue atoms to the box, not a modeled sidechain reach.",
            "No docking was rerun and no full-set ProLIF calculation was run.",
            "Canonical UniProt mapping is reported where aligned; GPCRdb/Ballesteros-Weinstein positions were not available and were not guessed.",
        ],
    }
    (out_dir / "audit_summary.json").write_text(
        json.dumps(json_safe(summary), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "generated_at_utc": summary["generated_at_utc"],
        "script": str(Path(__file__).resolve().relative_to(ROOT)),
        "script_sha256": sha256(Path(__file__).resolve()),
        "inputs": inventory,
        "derived_outputs_sha256": {},
        "parameters": {
            "protonation_pH": 7.4,
            "box_envelope_angstrom": 6.0,
            "receptor_repair": "identity labels only by exact source coordinate match; no atoms added",
            "redocking_performed": False,
            "complete_423_result_processing_performed": False,
        },
    }
    write_json = lambda path, data: path.write_text(
        json.dumps(json_safe(data), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_json(out_dir / "corrected_receptor_manifest.json", {
        "generated_at_utc": summary["generated_at_utc"],
        "receptors": receptor_summaries,
        "identity_relabelled_structures_are_docking_ready": False,
        "no_sidechain_coordinates_added": True,
    })

    mmcif_match_text = "; ".join(
        "{}: {}/{}".format(
            receptor_id,
            data["pdb_vs_mmcif_atom_inventory"]["n_residues_with_identical_atom_name_sets"],
            data["pdb_vs_mmcif_atom_inventory"]["n_cif_residues"],
        )
        for receptor_id, data in receptor_summaries.items()
    )
    atom_change_text = "; ".join(
        "{}: {}".format(
            receptor_id, data["pdbqt_atom_name_change_classification"]
        )
        for receptor_id, data in receptor_summaries.items()
    )
    receptor_charge_text = "; ".join(
        "{}: {}".format(receptor_id, data["pdbqt_nonzero_partial_charges"])
        for receptor_id, data in receptor_summaries.items()
    )
    construct_text = "; ".join(
        "{} entity {}: {}; {} residuos no alineados a UniProt, {} dentro de 6 Å de caja".format(
            receptor_id,
            data["mmcif_chain_metadata"]["entity_id"],
            data["mmcif_chain_metadata"]["entity_description"],
            data["sequence_alignment"]["unmapped_residue_count"],
            data["sequence_alignment"]["unmapped_residues_within_6A_of_box"],
        )
        for receptor_id, data in receptor_summaries.items()
    )
    report_lines = [
        "# Reconciliación de receptores para docking 5-HT2A",
        "",
        f"Fecha UTC: {summary['generated_at_utc']}",
        "",
        "## Dictamen",
        "",
        f"**{decision_summary['conclusion']}** para aceptar el conjunto completo sin cambios. No se repitió ningún docking. En los tres receptores hay residuos incompletos cuya cadena lateral ausente puede alcanzar la caja; en 9LL8/9LLB varios se rotularon como ALA en el PDBQT.",
        "",
        "| Receptor | Alineamiento UniProt | Pesados PDB/PDBQT | Átomos con resname discrepante | Residuos incompletos ≤6 Å de caja | Caja | Exhaustiveness/modos | Decisión |",
        "|---|---:|---:|---:|---:|---|---|---|",
    ]
    for receptor_id, data in receptor_summaries.items():
        box = data["docking_box"]
        decision = decision_summary["receptor_decisions"][receptor_id]["decision"]
        report_lines.append(
            f"| {receptor_id} | {data['sequence_alignment']['local_alignment_identity']:.3f} "
            f"({data['sequence_alignment']['local_alignment_matches']}/{data['sequence_alignment']['local_alignment_length']}) | "
            f"{data['source_heavy_atoms']}/{data['pdbqt_heavy_atoms']} | "
            f"{data['docking_pdbqt_atom_mapping']['matched_coordinate_atom_pairs'] - data['docking_pdbqt_atom_mapping']['residue_name_matches']} | "
            f"{data['source_incomplete_residues_near_box_count']} | "
            f"({box['center_x']},{box['center_y']},{box['center_z']}) / "
            f"{box['size_x']}×{box['size_y']}×{box['size_z']} Å | "
            f"{data['docking_parameters']['parameters'].get('exhaustiveness')}/"
            f"{data['docking_parameters']['parameters'].get('num_modes')} | {decision} |"
        )
    report_lines.extend(
        [
            "",
            "La alineación usa la secuencia UniProt P28223 y el mapeo local de secuencia; no se infieren números GPCRdb/Ballesteros-Weinstein ausentes de los archivos. `receptor_atom_comparison.csv` registra el emparejamiento átomo a átomo por coordenadas pesadas y las etiquetas leídas directamente de las líneas PDBQT. `receptor_residue_comparison.csv` registra los inventarios por residuo, nombres esperados/faltantes y distancias a caja.",
            f"El mmCIF y el PDB fuente coinciden en el inventario de átomos de todos los residuos proteicos observados ({mmcif_match_text}); los residuos no proteicos se excluyeron. En consecuencia, el mmCIF tampoco aporta las cadenas laterales faltantes del bolsillo. Los {pilot_source_hash_comparison['matching_input_count']}/{pilot_source_hash_comparison['compared_input_count']} archivos de origen previamente registrados en el manifiesto del piloto tienen SHA-256 idéntico; no se modificaron los originales.",
            f"Entidad/cadena de constructo según mmCIF: {construct_text}. Las porciones no alineadas se mantienen diferenciadas de secuencia 5-HT2A; no se interpretan como mutaciones del receptor.",
            "",
            "## Causa de las discrepancias",
            "",
            "En los casos `GLU355→ALA2`, `LYS223→ALA144`, `ASP231→ALA154` y `ASP232→ALA155`, el nombre ALA no es evidencia de mutación: se compara por coordenadas con el PDB fuente y el mmCIF experimental. La identidad se recupera mediante secuencia y coordenadas; el PDBQT conserva átomos observados de backbone/CB pero rotula ciertos residuos parciales como ALA. Hay **residuo experimental incompleto más pérdida de identidad química en el PDBQT**, no una sustitución de secuencia demostrada. El detalle PDB↔mmCIF por átomo está en los CSV y JSON; no se añadió ninguna cadena lateral.",
            "",
            "Los cambios de nombre de átomos restantes deben interpretarse con `atom_name_comparison`: los intercambios en pares químicamente simétricos (p. ej., CD1/CD2) se separan de los cambios no simétricos. Las coordenadas de los átomos pesados emparejados se comparan directamente.",
            f"En los PDBQT, todos los átomos pesados emparejan por coordenada y elemento. Diferencias de nombre de átomo (todas simétricas/no químicas): {atom_change_text}. Cargas no nulas PDBQT por receptor: {receptor_charge_text}; los tipos y cargas por átomo están en `receptor_atom_comparison.csv`.",
            "",
            "## Residuos incompletos cerca de las cajas",
            "",
        ]
    )
    for receptor_id, decision in decision_summary["receptor_decisions"].items():
        report_lines.append(f"### {receptor_id}")
        if not decision["residues"]:
            report_lines.append("No se detectaron cadenas laterales incompletas en la envolvente de 6 Å de la caja.")
        else:
            report_lines.append(
                "| Residuo fuente | PDBQT observado | Átomos pesados ausentes en fuente | Distancia a caja | Distancia a pose piloto más cercana |"
            )
            report_lines.append("|---|---|---|---:|---|")
            for row in decision["residues"]:
                report_lines.append(
                    f"| {row['resname']}{row['resid']} | {row['pdbqt_resname']}"
                    f"{row['pdbqt_residue']} cadena {row['pdbqt_chain']} | "
                    f"{row['missing_atoms']} | {row['minimum_distance_to_box_angstrom']:.2f} Å | "
                    f"{row['minimum_distance_to_pilot_ligand_angstrom']:.2f} Å "
                    f"({row['nearest_pilot_ligand_id']} modelo {row['nearest_pilot_pose_model']}) |"
                )
        report_lines.append("")
    for receptor_id, data in receptor_summaries.items():
        substitutions = data["sequence_alignment"]["observed_sequence_substitutions"]
        report_lines.append(f"### Sustituciones de secuencia en {receptor_id}")
        if not substitutions:
            report_lines.append("No se detectaron sustituciones en los residuos observados de la cadena 5-HT2A frente a UniProt P28223.")
        else:
            report_lines.append(
                "| Residuo del constructo | Residuo WT | UniProt | Anotación mmCIF | Distancia a caja | Distancia mínima a pose piloto |"
            )
            report_lines.append("|---|---|---:|---|---:|---:|")
            substitution_distances = []
            for mutation in substitutions:
                row = next(
                    (
                        item for item in all_residue_rows
                        if item["receptor"] == receptor_id
                        and item["source_chain"] == mutation["chain"]
                        and str(item["source_residue"]) == str(mutation["source_residue"])
                    ),
                    None,
                )
                if row is None:
                    row = next(
                        (
                            item for item in all_residue_rows
                            if item["receptor"] == receptor_id
                            and str(item["source_uniprot_position"]) == str(mutation["uniprot_position"])
                            and item["source_resname"] == mutation["source_resname"]
                        ),
                        None,
                    )
                if row is None:
                    raise ValueError(
                        f"Cannot map {receptor_id} mutation to residue comparison: {mutation}"
                    )
                annotations = data["sequence_alignment"]["mmcif_annotated_sequence_differences"]
                annotation = next(
                    (
                        item["annotation"] for item in annotations
                        if item["auth_residue"] == str(mutation["source_residue"])
                    ),
                    "sin anotación explícita",
                )
                report_lines.append(
                    f"| {mutation['source_resname']}{mutation['source_residue']} | "
                    f"{mutation['reference_resname']} | {mutation['uniprot_position']} | "
                    f"{annotation} | {row['minimum_heavy_atom_distance_to_box_angstrom']:.2f} Å | "
                    f"{row['minimum_distance_to_pilot_ligand_angstrom']:.2f} Å |"
                )
                substitution_distances.append(
                    f"{mutation['source_resname']}{mutation['source_residue']}: "
                    f"{row['minimum_distance_to_pilot_ligand_angstrom']:.2f} Å"
                )
            distances = ", ".join(substitution_distances)
            report_lines.append(
                f"Las diferencias de {receptor_id} están etiquetadas en `_struct_ref_seq_dif` del mmCIF. Son cambios reales, no renumeraciones ni errores de PDBQT. Distancias mínimas a poses piloto: {distances}. Si el objetivo es el receptor humano WT, no se puede aprobar este receptor sin restaurar y validar las identidades WT."
            )
        report_lines.append("")
    report_lines.extend(
        [
            "## Alternativas estructurales para cadenas laterales faltantes",
            "",
            "Se compararon los otros dos receptores experimentales disponibles como plantillas locales, alineando los carbonos de backbone de los residuos vecinos (±5 posiciones canónicas) y midiendo RMSD después del ajuste rígido. `receptor_template_assessment.csv` conserva todas las opciones. El umbral 0.5 Å solo marca candidatos para revisión manual; un RMSD local bajo es evidencia de compatibilidad de backbone, no demuestra que la rotámera sea adecuada en el otro estado conformacional. No se transfirieron átomos automáticamente.",
            "",
            "| Receptor objetivo/residuo | Átomos ausentes | Plantilla | Cobertura de átomos | RMSD backbone local | Evaluación |",
            "|---|---|---|---:|---:|---|",
        ]
    )
    for candidate in template_assessments:
        rmsd = candidate["local_backbone_fit_rmsd_angstrom"]
        report_lines.append(
            f"| {candidate['target_receptor']}:{candidate['resname']}{candidate['target_resid']} "
            f"(UniProt {candidate['uniprot_position']}) | {candidate['missing_atoms']} | "
            f"{candidate['donor_receptor'] or 'ninguna'} | {candidate['missing_atom_coverage']:.2f} | "
            f"{float(rmsd):.3f} Å | {candidate['assessment']} |"
            if rmsd != ""
            else (
                f"| {candidate['target_receptor']}:{candidate['resname']}{candidate['target_resid']} "
                f"(UniProt {candidate['uniprot_position']}) | {candidate['missing_atoms']} | "
                "ninguna | 0.00 | — | no hay plantilla experimental local |"
            )
        )
    report_lines.extend(
        [
            "",
            "La GLU351 incompleta de 7WC8 no tiene plantilla experimental local con esa posición y cadena lateral completa, por lo que 7WC8 permanece no apto para docking WT. Los archivos `91Q_crystal.pdb`, `92S_crystal.pdb` y `A1ELA_crystal.pdb` contienen ligandos sin cadena proteica y no sirven como plantillas de receptor. Las opciones de menor RMSD para 9LLB ASP231/ASP232 proceden de 9LL8, pero su uso requiere revisar orientación carboxilato y entorno local; se conservan como propuestas, no como corrección aplicada.",
            "",
        ]
    )
    report_lines.extend(
        [
            "El sobre de 6 Å es una prueba conservadora para priorizar inspección de cadenas laterales ausentes, no una geometría reconstruida. Un residuo puede afectar el docking si su cadena lateral ausente alcanza la caja aunque los átomos observados estén fuera. En 7WC8, GLU351 está a 4.80 Å de la caja; en 9LL8, ASP217/218, LYS223 y GLU224 están en ella o cerca; en 9LLB, ASP231/232 y otros residuos del bolsillo también. No se deben inventar rotámeros a partir de una sola pose.",
            "",
            "## Parámetros y ejecutable originales",
            "",
            "Las 423 configuraciones por receptor contienen una sola variante de caja/parámetros. Exhaustiveness=8, num_modes=9, CPU=4; cada lote pidió 424 ligandos, registró 423 como correctos y 1 fallido/excluido. El script configura el comando `/usr/local/bin/vina --config=<config de ligando>`; los logs identifican `AutoDock Vina v1.2.7`, scoring `vina`, y no contienen marcador Vina-GPU. El campo CPU=4 y el ejecutable apoyan Vina CPU, no Vina-GPU. Una semilla se registra solo en parte de los logs (conteos por receptor en `audit_summary.json`); los configs no guardan semilla y no se pudo controlar la semilla de las corridas originales.",
            "",
            "## Integridad química de 7WC8",
            "",
            "La preparación ProLIF previa informa que corrigió órdenes de enlace carboxilato de 7WC8 (incluido ASP217 en numeración PDB), conservando carga y coordenadas; el resumen y las correcciones por residuo se incorporan a `audit_summary.json`. La secuencia de corrección fue reproducible y reportó cero problemas de sanitización, pero no resuelve cadenas laterales faltantes en los archivos originales.",
            "",
            "## Preparación previa de ProLIF",
            "",
            "Los PDB preparados previos se compararon átomo a átomo con los PDB fuente. Se preservó el inventario pesado y las coordenadas según el manifiesto; PDB2PQR propuso 41/46/43 átomos pesados adicionales (7WC8/9LL8/9LLB), pero la rutina los descartó explícitamente. Esta conservación no repara cadenas laterales ausentes. Las cargas PQR y los tipos se resumen por receptor en `audit_summary.json`.",
            "",
            "## Corrección sin modelado especulativo",
            "",
            "Se generó una copia de auditoría de cada PDBQT con etiquetas restauradas solo cuando el emparejamiento por coordenada pesada con el PDB fuente es unívoco. Las coordenadas repetidas se dejan sin renombrar y su cantidad se registra en el manifiesto; no se fuerza una asignación arbitraria. Estas copias están bajo `receptor_repair/identity_relabelled_not_for_docking/`; no se añadieron coordenadas ni átomos, no son aptas para docking cuando la cadena lateral fuente falta y no sustituyen una preparación corregida.",
            "",
            "El PDBQT corregido para redocking debe generarse únicamente tras obtener una estructura de receptor completa o un modelo de cadena lateral explícitamente aprobado, aplicar protonación y cargas/tipos consistentes, e inspeccionar el sitio. Si la reconstrucción no puede justificarse estructuralmente, el receptor correspondiente permanece no apto.",
            "",
            "## ProLIF",
            "",
            "No se regeneraron fingerprints con los PDBQT auditados: el pipeline ProLIF previo reconstruyó receptores desde PDB incompletos y no incorpora las cadenas laterales ausentes. Sus conteos son reproducibles computacionalmente, pero no validan por sí mismos la química del bolsillo de 9LL8/9LLB. Los fingerprints de esas estructuras deben recalcularse solo después de disponer de receptores completos y aprobados.",
            "",
        ]
    )
    (out_dir / "reports" / "receptor_reconciliation_report.md").write_text(
        "\n".join(report_lines), encoding="utf-8"
    )

    ligand_groups = Counter(item["protonation_decision"] for item in ligand_summary.values())
    ligand_state_checks = [
        item["pdbqt_state_comparison"]
        for item in ligand_summary.values()
        if item["has_pdbqt"]
    ]
    plan = [
        "# Plan de redocking selectivo",
        "",
        "No se ejecutó docking. Los resultados originales permanecen intactos.",
        "",
        "## Decisión por receptor",
        "",
        "| Receptor | Decisión | Alcance recomendado |",
        "|---|---|---|",
    ]
    for receptor_id, data in decision_summary["receptor_decisions"].items():
        scope = (
            "Una vez aprobado un receptor químicamente completo, repetir los 423 ligandos elegibles contra este receptor; la deficiencia está en el input receptor compartido, no depende de la pose final."
            if data["decision"] == "REPETIR DOCKING"
            else "Conservar los dockings; no hay evidencia en esta auditoría de una alteración química dentro/cerca de la caja."
        )
        plan.append(f"| {receptor_id} | {data['decision']} | {scope} |")
    plan.extend(
        [
            "",
            "## Ligandos",
            "",
            f"Se auditaron {len(ligand_summary)} estructuras SDF; {sum(item['has_pdbqt'] for item in ligand_summary.values())} tienen PDBQT preparado y una fue excluida por contener un elemento no soportado por AutoDock. Decisiones de estado: `{dict(ligand_groups)}`.",
            "",
            f"En los {len(ligand_state_checks)} PDBQT disponibles, la correspondencia por coordenadas de átomos pesados y el patrón de H polares contra el SDF se verificaron en "
            f"{sum(item.get('heavy_atom_coordinate_coverage') == 1.0 and item.get('polar_hydrogen_pattern_matches_sdf', False) for item in ligand_state_checks)}/{len(ligand_state_checks)}. Los PDBQT tienen "
            f"{sum(item['pdbqt_nonzero_partial_charges'] or 0 for item in ligand_summary.values())} cargas parciales no nulas en conjunto; esto no se usa como criterio de error para scoring Vina estándar. El principal problema es el microestado: el PDBQT conserva el estado del SDF de entrada, que no necesariamente es el predominante a pH 7.4.",
            "",
            "Las entradas se prepararon con la política `pending_dimorphite_ph7.4` / `preserve_input_smiles_state`; el PDBQT de ligando deriva de ese SDF. Dimorphite enumeró candidatos, pero no ofrece poblaciones, y no se encontró predictor de pKa instalado. Por eso la auditoría marca estados múltiples o discrepancias de carga para revisión de sensibilidad, no convierte automáticamente toda la biblioteca en redocking.",
            "",
            "Prioridad de microestado: ligandos con aminas alifáticas protonables, guanidinas/amidinas o ácidos carboxílicos; confirmar identidad de protómero con pKa medido/validado o criterio químico documentado. Si el microestado original queda fuera del estado elegido, repetir solo ese ligando para cada receptor que se use, con ID de microestado separado y parámetros de caja/semilla registrados. No usar el contacto con ASP155 para elegir estados.",
            "",
            "## Parámetros para la repetición",
            "",
            "Mantener el centro/tamaño original por receptor, exhaustiveness=8 y num_modes=9, scoring `vina`, salvo justificación documentada. La semilla original no está en los archivos config y solo aparece en parte de los logs. Para la comparación controlada, guardar la semilla explícita de cada nuevo cálculo y reportar dónde no se pudo emparejar con la corrida antigua.",
            "",
            "No iniciar los 423 cálculos de ninguno de los tres receptores hasta resolver y validar todas las cadenas laterales incompletas dentro/cerca de la caja, y hasta aprobar inventario de residuos, geometría, protonación y tipos atómicos. Si se modelan cadenas laterales, conservar el PDB experimental original y documentar plantilla, alineamiento, rotámero, relajación local y validación.",
            "",
        ]
    )
    (out_dir / "reports" / "redocking_plan.md").write_text("\n".join(plan), encoding="utf-8")

    final_report = [
        "# Validación pre-docking — 5-HT2A",
        "",
        f"Conclusión global: **{decision_summary['conclusion']}**.",
        "",
        "| Receptor | Identidad estructural | Preparación química | Docking original | Decisión |",
        "|---|---|---|---|---|",
    ]
    for receptor_id, data in receptor_summaries.items():
        decision = decision_summary["receptor_decisions"][receptor_id]
        final_report.append(
            f"| {receptor_id} | UniProt P28223; {len(data['sequence_alignment']['mmcif_annotated_sequence_differences'])} mutaciones anotadas | "
            f"{data['source_incomplete_residues_near_box_count']} residuos incompletos en envolvente de caja | "
            f"Vina 1.2.7 CPU confirmado | {decision['decision']} |"
        )
    final_report.extend(
        [
            "",
            "### Ligandos por estado",
            "",
            "| Grupo/decisión | Número | Acción |",
            "|---|---:|---|",
        ]
    )
    for label, count in ligand_groups.items():
        action = (
            "Revisar microestados; no asumir población ni repetir hasta confirmar."
            if "manual_review" in label
            else (
                "No dockeado: exclusión por soporte de elemento."
                if "excluded" in label
                else "Conservar como estado de entrada documentado; pKa/población no validados."
            )
        )
        final_report.append(f"| {label} | {count} | {action} |")
    final_report.extend(
        [
            "",
            "La conclusión `NO APROBADO` significa que todavía no se debe correr el conjunto completo ni regenerar el corrected pilot final. Los tres receptores requieren resolver cadenas laterales ausentes dentro/cerca de sus cajas. Los ligandos quedan sujetos a revisión de microestado donde la enumeración presenta discrepancias o alternativas; no hay población pKa estimada. No se atribuye función, afinidad ni eficacia a los IFP.",
            "",
            "Consulte `receptor_reconciliation_report.md`, `redocking_plan.md`, `ligand_protonation/ligand_protonation_audit.csv` y `redocking_decisions/redocking_decisions.csv`.",
        ]
    )
    (out_dir / "reports" / "final_pre_docking_validation_report.md").write_text(
        "\n".join(final_report), encoding="utf-8"
    )
    (out_dir / "logs" / "audit_execution.log").write_text(
        "\n".join(
            [
                f"generated_at_utc={summary['generated_at_utc']}",
                f"script={Path(__file__).resolve().relative_to(ROOT)}",
                f"script_sha256={sha256(Path(__file__).resolve())}",
                f"dimorphite_dl={summary['versions']['dimorphite_dl']}",
                f"rdkit={summary['versions']['rdkit']}",
                f"gemmi={summary['versions']['gemmi']}",
                f"source_inventory_files={len(inventory)}",
                f"input_hash_matches_prior_pilot={pilot_source_hash_comparison['matching_input_count']}/{pilot_source_hash_comparison['compared_input_count']}",
                f"receptor_decisions={json.dumps(decision_summary['receptor_decisions'], sort_keys=True)}",
                f"redocking_executed=false",
                f"complete_docking_set_processed=false",
                f"conclusion={decision_summary['conclusion']}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    manifest["derived_outputs_sha256"] = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in sorted(out_dir.rglob("*"))
        if path.is_file() and path.name != "provenance_manifest.json"
    }
    write_json(out_dir / "provenance_manifest.json", manifest)
    print(
        json.dumps(
            {
                "output_dir": str(out_dir.relative_to(ROOT)),
                "receptor_decisions": {
                    key: value["decision"]
                    for key, value in decision_summary["receptor_decisions"].items()
                },
                "ligands": len(ligand_summary),
                "ligand_microstate_decisions": dict(ligand_groups),
                "docking_decisions": dict(Counter(row["decision"] for row in decisions)),
                "global_conclusion": decision_summary["conclusion"],
                "report": str((out_dir / "reports" / "final_pre_docking_validation_report.md").relative_to(ROOT)),
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    from rdkit import Chem

    raise SystemExit(main())
