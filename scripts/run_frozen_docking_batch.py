#!/usr/bin/env python3
"""Execute a preflighted docking manifest with input-hash checked resumption."""

from __future__ import annotations

import argparse
import csv
import csv
import hashlib
import json
import re
import subprocess
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


def refresh_manifest() -> tuple[list[dict[str, str]], dict[str, Any]]:
    summary_path = recovery.OUT / "recovery_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(
            "Recovery summary is missing; complete recovery and pilot validation first"
        )
    recovery_summary = read_json(summary_path)
    ligand_rows = recovery.read_csv(
        recovery.OUT / "ligand_preparation_summary.csv"
    )
    receptor_preparations = recovery_summary["receptor_preparations"]
    _, preflight = recovery.create_mass_docking_preflight(
        ligand_rows, receptor_preparations, recovery.OUT
    )
    manifest = recovery.read_csv(
        recovery.OUT / "mass_docking_preflight.csv"
    )
    return manifest, preflight


def run_manifest(manifest: list[dict[str, str]], preflight: dict[str, Any]) -> dict[str, int]:
    if preflight["status"] != "READY":
        raise RuntimeError(
            "Docking batch blocked: resolve all manual-review entries first; "
            f"status counts={preflight['status_counts']}"
        )
    executable = recovery.shutil.which("vina") or "/usr/local/sbin/vina"
    if not Path(executable).is_file():
        raise FileNotFoundError(f"AutoDock Vina executable not found: {executable}")
    version = subprocess.run(
        [executable, "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    if "1.2.7" not in version:
        raise RuntimeError(f"Expected AutoDock Vina 1.2.7, found {version!r}")

    ready_rows = [row for row in manifest if row["status"] == "READY"]
    batch_fingerprint = hashlib.sha256(
        json.dumps(
            [
                {
                    key: row[key]
                    for key in (
                        "ligand_id", "receptor", "ligand_sha256",
                        "receptor_sha256", "receptor_pdb_sha256", "box_sha256",
                        "ligand_summary_sha256", "recovery_script_sha256",
                        "prolif_script_sha256", "batch_script_sha256",
                        "ifp_script_sha256", "vina_executable_sha256",
                        "parameters_json", "seed",
                    )
                }
                for row in ready_rows
            ],
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()[:16]
    root = recovery.OUT / "final_docking" / batch_fingerprint
    counts = {"completed": 0, "resumed": 0, "failed": 0}
    result_rows = []
    for row in manifest:
        if row["status"] == "NO DISPONIBLE":
            result_rows.append({
                "receptor": row["receptor"],
                "ligand_id": row["ligand_id"],
                "status": "unavailable",
                "pose": "",
                "seed": row["seed"],
                "pose_count": 0,
                "top_affinity_kcal_mol": "",
                "error": row["reason"],
            })
    for row in ready_rows:
        ligand_path = ROOT / row["ligand_pdbqt"]
        receptor_path = ROOT / row["receptor_pdbqt"]
        box_path = ROOT / row["box_json"]
        if (
            file_hash(ligand_path) != row["ligand_sha256"]
            or file_hash(receptor_path) != row["receptor_sha256"]
            or file_hash(ROOT / row["receptor_pdb"])
            != row["receptor_pdb_sha256"]
            or file_hash(box_path) != row["box_sha256"]
            or file_hash(recovery.OUT / "ligand_preparation_summary.csv")
            != row["ligand_summary_sha256"]
            or file_hash(ROOT / "scripts" / "run_docking_recovery.py")
            != row["recovery_script_sha256"]
            or file_hash(ROOT / "scripts" / "run_prolif_corrected.py")
            != row["prolif_script_sha256"]
            or file_hash(ROOT / "scripts" / "run_frozen_docking_batch.py")
            != row["batch_script_sha256"]
            or file_hash(ROOT / "scripts" / "run_frozen_ifp.py")
            != row["ifp_script_sha256"]
            or file_hash(Path(executable)) != row["vina_executable_sha256"]
        ):
            raise RuntimeError(
                f"Frozen input hash mismatch for {row['receptor']}:{row['ligand_id']}"
            )
        parameters = json.loads(row["parameters_json"])
        box = read_json(box_path)
        run_dir = root / row["receptor"]
        pose_dir, log_dir, config_dir = (
            run_dir / name for name in ("poses", "logs", "configs")
        )
        for directory in (pose_dir, log_dir, config_dir):
            directory.mkdir(parents=True, exist_ok=True)
        pose_path = pose_dir / f"{row['ligand_id']}.pdbqt"
        meta_path = pose_dir / f"{row['ligand_id']}.meta.json"
        log_path = log_dir / f"{row['ligand_id']}.log"
        config_path = config_dir / f"{row['ligand_id']}.json"
        identity = {
            "ligand_sha256": row["ligand_sha256"],
            "receptor_sha256": row["receptor_sha256"],
            "receptor_pdb_sha256": row["receptor_pdb_sha256"],
            "box_sha256": row["box_sha256"],
            "ligand_summary_sha256": row["ligand_summary_sha256"],
            "recovery_script_sha256": row["recovery_script_sha256"],
            "prolif_script_sha256": row["prolif_script_sha256"],
            "batch_script_sha256": row["batch_script_sha256"],
            "ifp_script_sha256": row["ifp_script_sha256"],
            "vina_executable_sha256": row["vina_executable_sha256"],
            "parameters": parameters,
            "seed": int(row["seed"]),
            "vina_version": version,
        }
        if meta_path.is_file() and pose_path.is_file():
            metadata = read_json(meta_path)
            if (
                metadata.get("identity") == identity
                and metadata.get("pose_sha256") == file_hash(pose_path)
                and metadata.get("returncode") == 0
            ):
                counts["resumed"] += 1
                result_rows.append({
                    "receptor": row["receptor"],
                    "ligand_id": row["ligand_id"],
                    "status": "ok",
                    "pose": str(pose_path.relative_to(ROOT)),
                    "seed": row["seed"],
                    "pose_count": metadata.get("pose_count", 0),
                    "top_affinity_kcal_mol": metadata.get("top_affinity_kcal_mol", ""),
                    "error": "",
                })
                continue
            raise RuntimeError(
                f"Existing output has stale or unverifiable provenance: {pose_path}"
            )

        config = {
            **identity,
            "ligand": row["ligand_pdbqt"],
            "receptor": row["receptor_pdbqt"],
            "box": row["box_json"],
            "output": str(pose_path.relative_to(ROOT)),
        }
        config_path.write_text(
            json.dumps(config, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        temp_pose = pose_path.with_suffix(".pdbqt.tmp")
        command = [
            executable,
            "--receptor", str(receptor_path),
            "--ligand", str(ligand_path),
            "--out", str(temp_pose),
            "--center_x", str(box["center_x"]),
            "--center_y", str(box["center_y"]),
            "--center_z", str(box["center_z"]),
            "--size_x", str(box["size_x"]),
            "--size_y", str(box["size_y"]),
            "--size_z", str(box["size_z"]),
            "--exhaustiveness", str(parameters["exhaustiveness"]),
            "--num_modes", str(parameters["num_modes"]),
            "--energy_range", str(parameters["energy_range"]),
            "--cpu", str(parameters["cpu"]),
            "--seed", str(row["seed"]),
        ]
        result = subprocess.run(
            command, capture_output=True, text=True, timeout=1800
        )
        log_path.write_text(
            result.stdout + result.stderr, encoding="utf-8"
        )
        if result.returncode != 0 or not temp_pose.is_file():
            temp_pose.unlink(missing_ok=True)
            counts["failed"] += 1
            result_rows.append({
                "receptor": row["receptor"],
                "ligand_id": row["ligand_id"],
                "status": "failed",
                "pose": "",
                "seed": row["seed"],
                "pose_count": 0,
                "top_affinity_kcal_mol": "",
                "error": (result.stderr or result.stdout)[-500:],
            })
            continue
        score_match = re.search(
            r"^\s*1\s+(-?\d+\.\d+)", result.stdout, re.MULTILINE
        )
        if not score_match:
            temp_pose.unlink(missing_ok=True)
            counts["failed"] += 1
            result_rows.append({
                "receptor": row["receptor"],
                "ligand_id": row["ligand_id"],
                "status": "failed",
                "pose": "",
                "seed": row["seed"],
                "pose_count": 0,
                "top_affinity_kcal_mol": "",
                "error": "Vina returned no top affinity row",
            })
            continue
        temp_pose.replace(pose_path)
        pose_count = sum(
            line.startswith("MODEL ")
            for line in pose_path.read_text(encoding="utf-8").splitlines()
        ) or 1
        metadata = {
            "identity": identity,
            "returncode": result.returncode,
            "pose_sha256": file_hash(pose_path),
            "pose_count": pose_count,
            "top_affinity_kcal_mol": float(score_match.group(1)),
            "config": str(config_path.relative_to(ROOT)),
            "log": str(log_path.relative_to(ROOT)),
        }
        meta_path.write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        counts["completed"] += 1
        result_rows.append({
            "receptor": row["receptor"],
            "ligand_id": row["ligand_id"],
            "status": "ok",
            "pose": str(pose_path.relative_to(ROOT)),
            "seed": row["seed"],
            "pose_count": pose_count,
            "top_affinity_kcal_mol": float(score_match.group(1)),
            "error": "",
        })
    results_path = root / "final_docking_results.csv"
    results_path.parent.mkdir(parents=True, exist_ok=True)
    with results_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result_rows[0]))
        writer.writeheader()
        writer.writerows(result_rows)
    (root / "batch_status.json").write_text(
        json.dumps(
            {
                "batch_fingerprint": batch_fingerprint,
                "vina_version": version,
                "expected_runs": len(ready_rows),
                "counts": counts,
            },
            indent=2,
            ensure_ascii=False,
        ) + "\n",
        encoding="utf-8",
    )
    if counts["failed"]:
        raise RuntimeError(f"Docking failures recorded: {counts}")
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Run the complete ready batch; without this flag only refresh preflight.",
    )
    args = parser.parse_args()
    manifest, preflight = refresh_manifest()
    print(json.dumps(preflight, indent=2, ensure_ascii=False))
    if not args.execute:
        return 0
    counts = run_manifest(manifest, preflight)
    print(json.dumps(counts, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
