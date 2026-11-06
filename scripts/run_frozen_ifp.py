#!/usr/bin/env python3
"""Calculate ProLIF fingerprints for a completed, hash-verified docking batch."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_docking_recovery as recovery  # noqa: E402


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_prepared_ligands(
    ligand_rows: list[dict[str, str]],
    manifest_by_ligand: dict[str, dict[str, str]],
) -> dict[str, dict[str, Any]]:
    from run_prolif_corrected import add_protomer_coordinates, read_sdf_records, select_protomer

    records = read_sdf_records(recovery.LIGAND_SDF)
    prepared = {}
    for row in ligand_rows:
        ligand_id = row["ligand_id"]
        if not recovery.parse_bool(row["pdbqt_available"]):
            continue
        manifest_entry = manifest_by_ligand[ligand_id]
        pdbqt = ROOT / manifest_entry["ligand_pdbqt"]
        if file_hash(pdbqt) != manifest_entry["ligand_sha256"]:
            raise RuntimeError(f"Prepared ligand changed since docking: {ligand_id}")
        reference = records[ligand_id]
        protomer, mapping, report = recovery.prolif_pipeline.select_protomer(
            reference, recovery.PH
        )
        mol_3d = add_protomer_coordinates(reference, protomer, mapping)
        atom_map, indices, map_report = recovery.map_pdbqt_to_protomer(
            reference,
            mol_3d,
            mapping,
            pdbqt,
        )
        prepared[ligand_id] = {
            "reference": reference,
            "protomer": protomer,
            "protomer_mapping": mapping,
            "pdbqt": pdbqt,
            "converted": protomer,
            "atom_mapping": atom_map,
            "map_report": map_report,
            "prepared_indices": indices,
            "protonation_report": report,
        }
    return prepared


def execute(batch_fingerprint: str) -> dict[str, Any]:
    batch_root = recovery.OUT / "final_docking" / batch_fingerprint
    status_path = batch_root / "batch_status.json"
    results_path = batch_root / "final_docking_results.csv"
    if not status_path.is_file() or not results_path.is_file():
        raise FileNotFoundError("The frozen docking batch has not completed")
    batch_status = read_json(status_path)
    if batch_status["counts"]["failed"]:
        raise RuntimeError("The frozen docking batch contains failed runs")
    if (
        batch_status["counts"]["completed"] + batch_status["counts"]["resumed"]
        != batch_status["expected_runs"]
    ):
        raise RuntimeError("The frozen docking batch is incomplete")

    manifest = recovery.read_csv(recovery.OUT / "mass_docking_preflight.csv")
    ready_entry = next(row for row in manifest if row["status"] == "READY")
    frozen_scripts = (
        ("scripts/run_docking_recovery.py", "recovery_script_sha256"),
        ("scripts/run_prolif_corrected.py", "prolif_script_sha256"),
        ("scripts/run_frozen_docking_batch.py", "batch_script_sha256"),
        ("scripts/run_frozen_ifp.py", "ifp_script_sha256"),
    )
    for relative_path, hash_field in frozen_scripts:
        if file_hash(ROOT / relative_path) != ready_entry[hash_field]:
            raise RuntimeError(f"Pipeline script changed since docking: {relative_path}")
    if file_hash(recovery.OUT / "ligand_preparation_summary.csv") != ready_entry[
        "ligand_summary_sha256"
    ]:
        raise RuntimeError("Ligand preparation summary changed since docking")
    rows_by_pair = {
        (row["receptor"], row["ligand_id"]): row for row in manifest
    }
    docking_rows = recovery.read_csv(results_path)
    for row in docking_rows:
        if row["status"] != "ok":
            continue
        manifest_row = rows_by_pair[(row["receptor"], row["ligand_id"])]
        if (
            file_hash(ROOT / manifest_row["ligand_pdbqt"])
            != manifest_row["ligand_sha256"]
            or file_hash(ROOT / manifest_row["receptor_pdbqt"])
            != manifest_row["receptor_sha256"]
            or file_hash(ROOT / manifest_row["box_json"])
            != manifest_row["box_sha256"]
        ):
            raise RuntimeError(
                f"Docking inputs changed since execution: {row['receptor']}:{row['ligand_id']}"
            )
        pose_path = ROOT / row["pose"]
        meta_path = pose_path.with_suffix(".meta.json")
        metadata = read_json(meta_path)
        if metadata["pose_sha256"] != file_hash(pose_path):
            raise RuntimeError(f"Docking pose hash mismatch: {pose_path}")
    ligand_rows = recovery.read_csv(
        recovery.OUT / "ligand_preparation_summary.csv"
    )
    ligand_manifest = {
        row["ligand_id"]: row for row in manifest
        if row["status"] == "READY"
    }
    prepared = load_prepared_ligands(ligand_rows, ligand_manifest)
    pilot = [
        {
            "ligand_id": ligand_id,
            "label": next(
                row["pharmacological_label"]
                for row in ligand_rows if row["ligand_id"] == ligand_id
            ),
        }
        for ligand_id in prepared
    ]

    final_root = recovery.OUT / "final_ifp" / batch_fingerprint
    final_root.mkdir(parents=True, exist_ok=True)
    receptor_link = final_root / "receptor_preparation"
    if not receptor_link.exists():
        receptor_link.symlink_to(
            recovery.OUT / "receptor_preparation",
            target_is_directory=True,
        )
    pose_records, summary = recovery.process_pilot_prolif(
        pilot, prepared, docking_rows, final_root
    )
    result = {
        "batch_fingerprint": batch_fingerprint,
        "dockings": sum(row["status"] == "ok" for row in docking_rows),
        "prolif_poses": len(pose_records),
        "reproducibility": summary["reproducibility"],
        "receptors": summary["receptors"],
        "interaction_categories": summary["interaction_categories"],
    }
    (final_root / "final_ifp_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-fingerprint", required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if not args.execute:
        raise SystemExit("Pass --execute to calculate the full ProLIF batch")
    print(json.dumps(execute(args.batch_fingerprint), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
