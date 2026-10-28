"""Generate and validate one 3D SDF conformer per eligible docking ligand."""
import csv
import datetime as dt
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

try:
    from rawstore import sha256_file
except ModuleNotFoundError:
    from scripts.rawstore import sha256_file

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "processed" / "docking_ligands.csv"
STRUCTURES = ROOT / "structures"
OUTPUT = STRUCTURES / "docking_ligands_3d.sdf"
QC_OUTPUT = STRUCTURES / "docking_ligands_sdf_qc.csv"
MANIFEST = ROOT / "results" / "dataset_manifest.json"
RANDOM_SEED = 20261028
MAX_MINIMIZATION_ITERATIONS = 1000
SDF_PROPERTIES = (
    "ligand_id", "inchikey", "canonical_smiles", "original_smiles",
    "pharmacological_label", "pharmacological_subtype", "source_identifiers",
    "evidence_record_ids", "reference_ids", "protonation_policy",
    "conformer_method", "minimization_method", "minimization_energy",
    "minimization_status", "minimization_max_iterations", "random_seed",
)


def write_csv(path, rows):
    columns = list(rows[0]) if rows else [
        "ligand_id", "status", "reason", "minimization_method",
        "minimization_energy", "minimization_status",
        "minimization_max_iterations", "random_seed",
    ]
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def get_force_field(mol):
    if AllChem.MMFFHasAllMoleculeParams(mol):
        properties = AllChem.MMFFGetMoleculeProperties(
            mol, mmffVariant="MMFF94s"
        )
        if properties is not None:
            force_field = AllChem.MMFFGetMoleculeForceField(mol, properties)
            if force_field is not None:
                return force_field, "MMFF94s"
    if AllChem.UFFHasAllMoleculeParams(mol):
        force_field = AllChem.UFFGetMoleculeForceField(mol)
        if force_field is not None:
            return force_field, "UFF"
    raise ValueError("neither MMFF94s nor UFF has complete parameters")


def identity_from_mol(mol):
    no_h = Chem.RemoveHs(Chem.Mol(mol))
    smiles = Chem.MolToSmiles(no_h, canonical=True, isomericSmiles=True)
    return smiles, Chem.MolToInchiKey(no_h)


def validate_coordinates(mol):
    if mol.GetNumConformers() != 1:
        raise ValueError(f"expected exactly one conformer, found {mol.GetNumConformers()}")
    conformer = mol.GetConformer()
    for atom_index in range(mol.GetNumAtoms()):
        point = conformer.GetAtomPosition(atom_index)
        if not all(math.isfinite(value) for value in (point.x, point.y, point.z)):
            raise ValueError(f"non-finite coordinate at atom {atom_index}")


def prepare_one(row):
    ligand_id = row["ligand_id"]
    mol = Chem.MolFromSmiles(row["canonical_smiles"])
    if mol is None:
        raise ValueError("canonical_smiles cannot be parsed")
    smiles, inchikey = identity_from_mol(mol)
    if smiles != row["canonical_smiles"]:
        raise ValueError("canonical SMILES identity mismatch")
    if inchikey != row["inchikey"] or f"LIG_{inchikey}" != ligand_id:
        raise ValueError("InChIKey/ligand_id mismatch")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("multiple fragments are outside the approved SDF set")

    original = Chem.MolFromSmiles(row["original_smiles"])
    if original is None:
        raise ValueError("original SMILES cannot be parsed")
    _, original_inchikey = identity_from_mol(original)
    if original_inchikey != inchikey:
        raise ValueError("original SMILES does not match canonical ligand identity")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    params.numThreads = 1
    conformer_id = AllChem.EmbedMolecule(mol, params)
    if conformer_id < 0:
        raise ValueError("ETKDGv3 failed to generate a conformer")
    validate_coordinates(mol)

    force_field, method = get_force_field(mol)
    iterations = force_field.Minimize(maxIts=MAX_MINIMIZATION_ITERATIONS)
    energy = float(force_field.CalcEnergy())
    if iterations != 0:
        raise ValueError(
            f"{method} minimization did not converge within "
            f"{MAX_MINIMIZATION_ITERATIONS} iterations (status={iterations})"
        )
    if not math.isfinite(energy):
        raise ValueError("minimization returned a non-finite energy")
    validate_coordinates(mol)
    output_smiles, output_inchikey = identity_from_mol(mol)
    if output_smiles != row["canonical_smiles"] or output_inchikey != inchikey:
        raise ValueError("3D generation changed the molecular identity")

    properties = {
        "ligand_id": ligand_id,
        "inchikey": inchikey,
        "canonical_smiles": row["canonical_smiles"],
        "original_smiles": row["original_smiles"],
        "pharmacological_label": row["pharmacological_label"],
        "pharmacological_subtype": row["pharmacological_subtype"],
        "source_identifiers": row["source_identifiers"],
        "evidence_record_ids": row["evidence_record_ids"],
        "reference_ids": row["reference_ids"],
        "protonation_policy": "preserve_input_smiles_state",
        "conformer_method": "RDKit_ETKDGv3",
        "minimization_method": method,
        "minimization_energy": f"{energy:.8f}",
        "minimization_status": "converged",
        "minimization_max_iterations": str(MAX_MINIMIZATION_ITERATIONS),
        "random_seed": str(RANDOM_SEED),
    }
    for name in SDF_PROPERTIES:
        mol.SetProp(name, properties[name])
    return mol, {
        "ligand_id": ligand_id,
        "status": "validated",
        "reason": "",
        "minimization_method": method,
        "minimization_energy": f"{energy:.8f}",
        "minimization_status": "converged",
        "minimization_max_iterations": MAX_MINIMIZATION_ITERATIONS,
        "random_seed": RANDOM_SEED,
    }


def validate_written_sdf(path, expected_rows):
    seen = set()
    errors = []
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    for index, mol in enumerate(supplier, start=1):
        if mol is None:
            errors.append({"record": index, "reason": "RDKit could not parse SDF record"})
            continue
        ligand_id = mol.GetProp("ligand_id") if mol.HasProp("ligand_id") else ""
        if not ligand_id or ligand_id in seen:
            errors.append({"record": index, "reason": f"missing/duplicate ligand_id:{ligand_id}"})
            continue
        seen.add(ligand_id)
        expected = expected_rows.get(ligand_id)
        if expected is None:
            errors.append({"record": index, "ligand_id": ligand_id, "reason": "unexpected ligand"})
            continue
        try:
            smiles, inchikey = identity_from_mol(mol)
            if smiles != expected["canonical_smiles"]:
                raise ValueError("round-trip canonical SMILES mismatch")
            if inchikey != expected["inchikey"]:
                raise ValueError("round-trip InChIKey mismatch")
            if mol.GetNumAtoms() <= Chem.RemoveHs(mol).GetNumAtoms():
                raise ValueError("no explicit hydrogens in SDF")
            validate_coordinates(mol)
            for property_name in SDF_PROPERTIES:
                if not mol.HasProp(property_name):
                    raise ValueError(f"missing SDF property {property_name}")
        except Exception as exc:
            errors.append({"record": index, "ligand_id": ligand_id, "reason": str(exc)})
    missing = set(expected_rows) - seen
    errors.extend(
        {"ligand_id": ligand_id, "reason": "missing from SDF output"}
        for ligand_id in sorted(missing)
    )
    if len(seen) != len(expected_rows):
        errors.append({
            "reason": f"record count mismatch: expected {len(expected_rows)}, found {len(seen)}"
        })
    return errors


def update_manifest(files, input_path):
    if MANIFEST.exists():
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    else:
        manifest = {}
    entries = {row["path"]: row for row in manifest.get("files", [])}
    for path in files:
        rel = path.relative_to(ROOT).as_posix()
        row_count = None
        if path.suffix.lower() == ".csv":
            with path.open(encoding="utf-8", newline="") as stream:
                row_count = max(sum(1 for _ in csv.reader(stream)) - 1, 0)
        elif path.suffix.lower() == ".sdf":
            with path.open(encoding="utf-8", errors="replace") as stream:
                row_count = sum(line.strip() == "$$$$" for line in stream)
        entries[rel] = {
            "path": rel,
            "rows": row_count,
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    inputs = manifest.get("inputs", [])
    input_rel = input_path.relative_to(ROOT).as_posix()
    input_record = {"path": input_rel, "sha256": sha256_file(input_path)}
    by_input = {item["path"]: item for item in inputs}
    by_input[input_rel] = input_record
    manifest["generated_at"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    manifest["versions"] = {
        **manifest.get("versions", {}),
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "rdkit": rdBase.rdkitVersion,
    }
    manifest["inputs"] = sorted(by_input.values(), key=lambda item: item["path"])
    manifest["files"] = sorted(entries.values(), key=lambda item: item["path"])
    note = (
        "Ligand SDF generated with RDKit ETKDGv3 and minimized with MMFF94s "
        "(UFF fallback), preserving input SMILES protonation state."
    )
    notes = manifest.setdefault("notes", [])
    if note not in notes:
        notes.append(note)
    temporary = MANIFEST.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(MANIFEST)


def main():
    if not INPUT.exists():
        print(f"ERROR: missing docking ligand input: {INPUT}", file=sys.stderr)
        return 2
    if OUTPUT.exists() or QC_OUTPUT.exists():
        print(
            f"ERROR: output exists; refusing to overwrite: {OUTPUT}, {QC_OUTPUT}",
            file=sys.stderr,
        )
        return 3
    with INPUT.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        print("ERROR: docking ligand CSV is empty", file=sys.stderr)
        return 4
    if any(row["docking_eligibility"] != "eligible_for_structure_preparation" for row in rows):
        print("ERROR: input CSV contains a molecule not cleared for preparation", file=sys.stderr)
        return 5

    STRUCTURES.mkdir(parents=True, exist_ok=True)
    output_tmp = OUTPUT.with_suffix(".sdf.tmp")
    writer = Chem.SDWriter(str(output_tmp))
    writer.SetForceV3000(False)
    qc_rows = []
    expected = {}
    for number, row in enumerate(rows, start=1):
        try:
            mol, qc = prepare_one(row)
            writer.write(mol)
            expected[row["ligand_id"]] = row
            qc_rows.append(qc)
        except Exception as exc:
            qc_rows.append({
                "ligand_id": row.get("ligand_id", ""),
                "status": "failed",
                "reason": f"{type(exc).__name__}: {exc}",
                "minimization_method": "",
                "minimization_energy": "",
                "minimization_status": "failed",
                "minimization_max_iterations": MAX_MINIMIZATION_ITERATIONS,
                "random_seed": RANDOM_SEED,
            })
        if number % 500 == 0:
            print(f"processed {number}/{len(rows)}", flush=True)
    writer.close()

    if any(row["status"] != "validated" for row in qc_rows):
        output_tmp.unlink(missing_ok=True)
        write_csv(QC_OUTPUT, qc_rows)
        print(
            f"ERROR: {sum(row['status'] != 'validated' for row in qc_rows)} "
            "ligands failed preparation; no SDF promoted to final path",
            file=sys.stderr,
        )
        return 6

    validation_errors = validate_written_sdf(output_tmp, expected)
    if validation_errors:
        output_tmp.unlink(missing_ok=True)
        qc_rows.append({
            "ligand_id": "",
            "status": "round_trip_validation_failed",
            "reason": json.dumps(validation_errors[:20], ensure_ascii=False),
            "minimization_method": "",
            "minimization_energy": "",
            "minimization_status": "round_trip_validation_failed",
            "minimization_max_iterations": MAX_MINIMIZATION_ITERATIONS,
            "random_seed": RANDOM_SEED,
        })
        write_csv(QC_OUTPUT, qc_rows)
        print(
            f"ERROR: SDF round-trip validation failed ({len(validation_errors)} issues)",
            file=sys.stderr,
        )
        return 7

    output_tmp.replace(OUTPUT)
    write_csv(QC_OUTPUT, qc_rows)
    update_manifest([OUTPUT, QC_OUTPUT], INPUT)
    methods = {}
    for row in qc_rows:
        methods[row["minimization_method"]] = methods.get(row["minimization_method"], 0) + 1
    print(json.dumps({
        "input_ligands": len(rows),
        "validated_sdf_records": len(expected),
        "failures": 0,
        "minimization_methods": methods,
        "sdf": OUTPUT.relative_to(ROOT).as_posix(),
        "qc": QC_OUTPUT.relative_to(ROOT).as_posix(),
        "seed": RANDOM_SEED,
        "protonation_policy": "preserve_input_smiles_state",
    }, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
