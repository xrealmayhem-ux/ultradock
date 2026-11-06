# 5-HT2A recovery, docking and ProLIF

This workflow is separate from the original UltraDock ligand-curation bundle. It preserves those inputs and records a frozen, hash-checked recovery run for 7WC8, 9LL8 and 9LLB.

## Frozen scientific inputs

- 424 SDF records were prepared from the original structures; 423 have usable original docking PDBQT inputs. One ligand is explicitly unavailable.
- A representative pH 7.4 state was selected for each of the 72 ambiguous ligands and approved by the project owner. Alternative states remain recorded in `results/docking_recovery/ligand_preparation/microstate_review_queue.csv`; populations were not estimated.
- 9LL8 LYS146 uses the 9LLB experimental side-chain template and ASP218 uses 7WC8. This paired selection removes the short contact seen when both were grafted from 7WC8. The corrected receptor and its hashes are included.
- The historical validation pilot used AutoDock Vina 1.2.7 on CPU. The final batch uses the existing Vina-GPU installation (`~/Vina-GPU` by default); set `ULTRADOCK_VINA_GPU_DIR` if it is elsewhere. The manifest freezes the GPU executable, wrapper and kernel hashes.
- Final-batch parameters: Vina scoring, 1000 GPU lanes, automatic search depth, 9 modes and energy range 3 kcal/mol. These settings are not equivalent to the CPU pilot's exhaustiveness 8. Per-combination seeds and box dimensions are in `results/docking_recovery/mass_docking_preflight.csv`.
- The 36 pilot dockings are historical and are not reused: they were run before the final 9LL8 correction or without per-run input hashes.

## Transferred checkpoint

`results/docking_recovery/final_docking/1ebb27eb0a58e56c/` contains 92 completed historical CPU dockings (pose, log, config and verified metadata). The full selection has 1,269 available ligand-receptor combinations; three combinations are explicitly unavailable. Those CPU poses are not reused by the new GPU batch.

## Verify the GPU runtime

The final runner requires the Vina-GPU executable, `run_vina.sh`, and `Kernel2_Opt.bin` in `~/Vina-GPU` (or `ULTRADOCK_VINA_GPU_DIR`). It also requires a working OpenCL GPU runtime. The preflight records hashes for these runtime files and will stop if they change. Verify the transferred prepared inputs with:

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

The runner uses Vina-GPU. The historical 92 CPU poses are not reused; the GPU batch gets a new fingerprint and runs all available combinations (1,269). Outputs are written under that fingerprint. A docking is resumable only when its input hashes, parameters, seed and Vina-GPU runtime match; incomplete temporary poses are not treated as completed results.

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
