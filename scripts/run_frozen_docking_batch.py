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


def vina_gpu_receptor_content(content: bytes) -> bytes:
    return b"".join(
        line
        for line in content.splitlines(keepends=True)
        if line.strip() != b"END"
    )


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
    ready_rows = [row for row in manifest if row["status"] == "READY"]
    if not ready_rows:
        raise RuntimeError("Docking batch has no ready combinations")
    engine_fields = (
        "vina_executable", "vina_executable_sha256",
        "vina_wrapper", "vina_wrapper_sha256",
        "vina_kernel", "vina_kernel_sha256", "vina_version",
    )
    engine_identity = {key: ready_rows[0][key] for key in engine_fields}
    if any(
        {key: row[key] for key in engine_fields} != engine_identity
        for row in ready_rows[1:]
    ):
        raise RuntimeError("Manifest contains inconsistent Vina-GPU runtimes")
    executable = engine_identity["vina_executable"]
    wrapper = engine_identity["vina_wrapper"]
    version = subprocess.run(
        [executable, "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    if version != engine_identity["vina_version"] or "Vina-GPU" not in version:
        raise RuntimeError(
            f"Vina-GPU version changed since preflight: "
            f"{engine_identity['vina_version']!r} -> {version!r}"
        )

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
                        "vina_wrapper_sha256", "vina_kernel_sha256",
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
            or file_hash(Path(wrapper)) != row["vina_wrapper_sha256"]
            or file_hash(Path(row["vina_kernel"])) != row["vina_kernel_sha256"]
        ):
            raise RuntimeError(
                f"Frozen input hash mismatch for {row['receptor']}:{row['ligand_id']}"
            )
        parameters = json.loads(row["parameters_json"])
        box = read_json(box_path)
        run_dir = root / row["receptor"]
        pose_dir, log_dir, config_dir, input_dir = (
            run_dir / name for name in ("poses", "logs", "configs", "inputs")
        )
        for directory in (pose_dir, log_dir, config_dir, input_dir):
            directory.mkdir(parents=True, exist_ok=True)
        gpu_receptor_path = input_dir / f"{row['receptor']}_vina_gpu.pdbqt"
        receptor_content = receptor_path.read_bytes()
        gpu_receptor_content = vina_gpu_receptor_content(receptor_content)
        if gpu_receptor_path.exists():
            if gpu_receptor_path.read_bytes() != gpu_receptor_content:
                raise RuntimeError(
                    f"Existing Vina-GPU receptor differs from frozen input: "
                    f"{gpu_receptor_path}"
                )
        else:
            gpu_receptor_path.write_bytes(gpu_receptor_content)
        gpu_receptor_sha256 = file_hash(gpu_receptor_path)
        pose_path = pose_dir / f"{row['ligand_id']}.pdbqt"
        meta_path = pose_dir / f"{row['ligand_id']}.meta.json"
        log_path = log_dir / f"{row['ligand_id']}.log"
        config_path = config_dir / f"{row['ligand_id']}.json"
        vina_config_path = config_dir / f"{row['ligand_id']}.txt"
        identity = {
            "ligand_sha256": row["ligand_sha256"],
            "receptor_sha256": row["receptor_sha256"],
            "gpu_receptor_sha256": gpu_receptor_sha256,
            "receptor_pdb_sha256": row["receptor_pdb_sha256"],
            "box_sha256": row["box_sha256"],
            "ligand_summary_sha256": row["ligand_summary_sha256"],
            "recovery_script_sha256": row["recovery_script_sha256"],
            "prolif_script_sha256": row["prolif_script_sha256"],
            "batch_script_sha256": row["batch_script_sha256"],
            "ifp_script_sha256": row["ifp_script_sha256"],
            "vina_executable_sha256": row["vina_executable_sha256"],
            "vina_wrapper_sha256": row["vina_wrapper_sha256"],
            "vina_kernel_sha256": row["vina_kernel_sha256"],
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
            "gpu_receptor": str(gpu_receptor_path.relative_to(ROOT)),
            "box": row["box_json"],
            "output": str(pose_path.relative_to(ROOT)),
        }
        config_path.write_text(
            json.dumps(config, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        temp_pose = pose_path.with_suffix(".pdbqt.tmp")
        vina_config_path.write_text(
            "\n".join(
                [
                    f"receptor = {gpu_receptor_path}",
                    f"ligand = {ligand_path}",
                    f"out = {temp_pose}",
                    f"center_x = {box['center_x']}",
                    f"center_y = {box['center_y']}",
                    f"center_z = {box['center_z']}",
                    f"size_x = {box['size_x']}",
                    f"size_y = {box['size_y']}",
                    f"size_z = {box['size_z']}",
                    f"thread = {parameters['thread']}",
                    f"search_depth = {parameters['search_depth']}",
                    f"num_modes = {parameters['num_modes']}",
                    f"energy_range = {parameters['energy_range']}",
                    f"seed = {row['seed']}",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        command = [wrapper, f"--config={vina_config_path}"]
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
                "engine": "Vina-GPU",
                "vina_version": version,
                "vina_executable_sha256": engine_identity["vina_executable_sha256"],
                "vina_kernel_sha256": engine_identity["vina_kernel_sha256"],
                "parameters": json.loads(ready_rows[0]["parameters_json"]),
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
