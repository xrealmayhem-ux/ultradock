import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import run_docking_recovery as recovery
import run_frozen_docking_batch as frozen_batch


class DockingRecoveryTests(unittest.TestCase):
    def test_vina_gpu_receptor_content_removes_end_marker_only(self):
        source = b"REMARK keep\r\nATOM atom\r\nEND\r\n"

        self.assertEqual(
            frozen_batch.vina_gpu_receptor_content(source),
            b"REMARK keep\r\nATOM atom\r\n",
        )

    def test_parse_bool_reads_csv_false_as_false(self):
        self.assertFalse(recovery.parse_bool("False"))
        self.assertFalse(recovery.parse_bool(""))
        self.assertTrue(recovery.parse_bool("True"))

    def test_create_decisions_repeats_changed_receptor_for_unambiguous_ligand(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            docking_dir = root / "original_docking"
            pose_dir = docking_dir / "poses"
            pose_dir.mkdir(parents=True)
            (pose_dir / "LIG-test.pdbqt").write_text("pose\n", encoding="utf-8")
            receptor = {
                "docking_dir": docking_dir,
            }
            ligand = {
                "ligand_id": "LIG-test",
                "preparation_status": "prepared",
                "state_ambiguity": "False",
                "error": "",
                "original_docking_state_differs_by_charge": "False",
                "original_formal_charge": 0,
                "selected_formal_charge": 0,
                "selected_state_smiles": "C",
            }
            with patch.dict(recovery.RECEPTORS, {"7WC8": receptor}, clear=True):
                decisions = recovery.create_decisions(
                    [ligand],
                    {"7WC8": {"repairs": [{"status": "grafted_from_template"}]}},
                    root,
                )
            self.assertEqual(decisions[0]["decision"], "REPETIR")

    def test_create_decisions_keeps_ambiguous_ligand_for_review(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            docking_dir = root / "original_docking"
            pose_dir = docking_dir / "poses"
            pose_dir.mkdir(parents=True)
            (pose_dir / "LIG-test.pdbqt").write_text("pose\n", encoding="utf-8")
            with patch.dict(
                recovery.RECEPTORS,
                {"7WC8": {"docking_dir": docking_dir}},
                clear=True,
            ):
                decisions = recovery.create_decisions(
                    [{
                        "ligand_id": "LIG-test",
                        "preparation_status": "prepared_state_requires_manual_review",
                        "state_ambiguity": "True",
                        "error": "",
                        "original_docking_state_differs_by_charge": "True",
                        "original_formal_charge": 0,
                        "selected_formal_charge": 1,
                        "selected_state_smiles": "[NH3+]C",
                    }],
                    {"7WC8": {"repairs": [{"status": "grafted_from_template"}]}},
                    root,
                )
            self.assertEqual(decisions[0]["decision"], "REVISIÓN MANUAL")


if __name__ == "__main__":
    unittest.main()
