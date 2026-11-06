# 5-HT2A recovery, docking and ProLIF

This workflow is separate from the original UltraDock ligand-curation bundle. It preserves those inputs and records a frozen, hash-checked recovery run for 7WC8, 9LL8 and 9LLB.

## Frozen scientific inputs

- 424 SDF records were prepared from the original structures; 423 have usable original docking PDBQT inputs. One ligand is explicitly unavailable.
- A representative pH 7.4 state was selected for each of the 72 ambiguous ligands and approved by the project owner. Alternative states remain recorded in `results/docking_recovery/ligand_preparation/microstate_review_queue.csv`; populations were not estimated.
- 9LL8 LYS146 uses the 9LLB experimental side-chain template and ASP218 uses 7WC8. This paired selection removes the short contact seen when both were grafted from 7WC8. The corrected receptor and its hashes are included.
- Vina 1.2.7 is bundled as a statically linked Linux x86-64 executable at `tools/vina`, matching the executable hash recorded in the manifest.
- Docking parameters: Vina scoring, exhaustiveness 8, 9 modes, energy range 3 kcal/mol and 4 CPUs. Per-combination seeds and box dimensions are in `results/docking_recovery/mass_docking_preflight.csv`.
- The 36 pilot dockings are historical and are not reused: they were run before the final 9LL8 correction or without per-run input hashes.

## Transferred checkpoint

`results/docking_recovery/final_docking/1ebb27eb0a58e56c/` contains 92 completed final dockings (pose, log, config and verified metadata). The complete run has 1,269 available ligand-receptor combinations; three combinations are explicitly unavailable. Thus 1,177 dockings remain.

## Recreate the environment

On Linux x86-64 with Conda or Mamba:

```bash
mamba env create -f environment-docking.yml
conda activate ultradock-5ht2a
export PATH="$PWD/tools:$PATH"
sha256sum tools/vina
```

The executable checksum must be:

```text
f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644
```

The frozen runner requires this Vina build and exact hashes for inputs, scripts and configuration. Do not edit the prepared structures or scripts after resuming.

On another Linux x86-64 machine, clone this repository normally; the 92 checked poses and prepared inputs are committed here. Verify the transferred snapshot with:

```bash
python scripts/verify_recovery_bundle.py
```

## Resume docking

First inspect the preflight; this does not launch Vina:

```bash
python scripts/run_frozen_docking_batch.py
```

Expected status is `READY`, with 1,269 ready and 3 unavailable combinations. Then resume:

```bash
python scripts/run_frozen_docking_batch.py --execute
```

The runner validates the checkpoint metadata, skips the 92 complete combinations, and runs the remaining 1,177. Outputs are written under the same fingerprinted directory. A docking is resumable only when its input hashes, parameters, seed and Vina version match; incomplete temporary poses are not treated as completed results.

## Calculate final IFP

After the batch reports all available runs complete:

```bash
python scripts/run_frozen_ifp.py \
  --batch-fingerprint 1ebb27eb0a58e56c \
  --execute
```

Fingerprints and per-receptor matrices are written under `results/docking_recovery/final_ifp/1ebb27eb0a58e56c/`. The final IFP runner checks docking pose hashes and the frozen input/script hashes before processing.

## Audit trail

- `results/docking_recovery/recovery_summary.json`: recovery tools, inputs and status.
- `results/docking_recovery/mass_docking_preflight.csv`: one row per ligand-receptor combination, input hashes, decisions, parameters and seeds.
- `results/docking_recovery/ligand_preparation_summary.csv`: selected states and atom-map checks.
- `results/docking_recovery/receptor_preparation/receptor_repair_manifest.csv`: localized graft evidence.
- The 92 checkpoint metadata files contain each pose hash, input identity, Vina version, score and output paths.
- `results/docking_recovery/transfer_manifest.json` covers the transferred snapshot with file sizes and SHA-256 hashes.

The full 423-ligand dataset and all 1,269 dockings must not be interpreted as experimental binding affinity or function.
