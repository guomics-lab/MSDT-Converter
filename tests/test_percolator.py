import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.percolator import (
    DuplicatePsmIdError,
    enrich_parquet_with_percolator,
    load_percolator_psms,
)


PERCOLATOR_COLUMNS = [
    "PSMId",
    "score",
    "q-value",
    "posterior_error_prob",
    "peptide",
    "proteinIds",
]


class PercolatorPsmTests(unittest.TestCase):
    def test_target_and_decoy_become_canonical_psm_records(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            target_path = root / "run_percolator_target_psms.tsv"
            decoy_path = root / "run_percolator_decoy_psms.tsv"
            pd.DataFrame(
                [["run.101.101.2_1", 5.2, 0.001, 0.002, "K.AC[57.0215]DM[15.9949].R", "P1"]],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(target_path, sep="\t", index=False)
            pd.DataFrame(
                [["run.102.102.3_1", -1.2, 0.4, 0.8, "R.PEPTIDE.K", "rev_P2"]],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(decoy_path, sep="\t", index=False)

            psms = load_percolator_psms(target_path, decoy_path)

            self.assertEqual(
                psms["psm_id"].tolist(),
                ["101_2_AC[57.02]DM[15.99]", "102_3_PEPTIDE"],
            )
            self.assertEqual(psms["label"].tolist(), [1, 0])
            self.assertEqual(psms["score"].tolist(), [5.2, -1.2])
            self.assertEqual(psms["q-value"].tolist(), [0.001, 0.4])
            self.assertEqual(psms["PEP"].tolist(), [0.002, 0.8])

    def test_duplicate_psm_id_across_target_and_decoy_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            target_path = root / "target.tsv"
            decoy_path = root / "decoy.tsv"
            duplicate = [
                "run.101.101.2_1",
                5.2,
                0.001,
                0.002,
                "K.PEPTIDE.R",
                "P1",
            ]
            pd.DataFrame([duplicate], columns=PERCOLATOR_COLUMNS).to_csv(
                target_path, sep="\t", index=False
            )
            pd.DataFrame([duplicate], columns=PERCOLATOR_COLUMNS).to_csv(
                decoy_path, sep="\t", index=False
            )

            with self.assertRaisesRegex(DuplicatePsmIdError, "101_2_PEPTIDE"):
                load_percolator_psms(target_path, decoy_path)

    def test_existing_parquet_is_strictly_matched_and_enriched(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_path = root / "input.parquet"
            output_path = root / "output.parquet"
            target_path = root / "target.tsv"
            decoy_path = root / "decoy.tsv"
            pd.DataFrame(
                {
                    "scan": [101, 102],
                    "charge": [2, 3],
                    "precursor_sequence": [
                        "AC[57.02]DM[15.99]",
                        "PEPTIDE",
                    ],
                    "label": [1, 0],
                    "precursor_mz": [500.2, 600.3],
                }
            ).to_parquet(input_path, index=False)
            pd.DataFrame(
                [["run.101.101.2_1", 5.2, 0.001, 0.002, "K.AC[57.0215]DM[15.9949].R", "P1"]],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(target_path, sep="\t", index=False)
            pd.DataFrame(
                [["run.102.102.3_1", -1.2, 0.4, 0.8, "R.PEPTIDE.K", "rev_P2"]],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(decoy_path, sep="\t", index=False)

            report = enrich_parquet_with_percolator(
                input_path, target_path, decoy_path, output_path
            )

            result = pd.read_parquet(output_path)
            self.assertEqual(result["scan"].tolist(), [101, 102])
            self.assertNotIn("psm_id", result.columns)
            self.assertEqual(result["score"].tolist(), [5.2, -1.2])
            self.assertEqual(result["q-value"].tolist(), [0.001, 0.4])
            self.assertEqual(result["PEP"].tolist(), [0.002, 0.8])
            self.assertEqual(report.input_rows, 2)
            self.assertEqual(report.matched_rows, 2)
            self.assertEqual(report.unmatched_input_rows, 0)

    def test_unmatched_parquet_psm_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_path = root / "input.parquet"
            target_path = root / "target.tsv"
            decoy_path = root / "decoy.tsv"
            pd.DataFrame(
                {
                    "scan": [101, 999],
                    "charge": [2, 2],
                    "precursor_sequence": ["PEPTIDE", "UNMATCHED"],
                    "label": [1, 1],
                }
            ).to_parquet(input_path, index=False)
            pd.DataFrame(
                [["run.101.101.2_1", 5.2, 0.001, 0.002, "K.PEPTIDE.R", "P1"]],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(target_path, sep="\t", index=False)
            pd.DataFrame(columns=PERCOLATOR_COLUMNS).to_csv(
                decoy_path, sep="\t", index=False
            )

            with self.assertRaisesRegex(
                ValueError, "MSDT Parquet PSMs missing from Percolator TSVs.*999_2_UNMATCHED"
            ):
                enrich_parquet_with_percolator(
                    input_path, target_path, decoy_path, root / "output.parquet"
                )

    def test_empty_percolator_metric_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            target_path = root / "target.tsv"
            decoy_path = root / "decoy.tsv"
            pd.DataFrame(
                [["run.101.101.2_1", 5.2, None, 0.002, "K.PEPTIDE.R", "P1"]],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(target_path, sep="\t", index=False)
            pd.DataFrame(columns=PERCOLATOR_COLUMNS).to_csv(
                decoy_path, sep="\t", index=False
            )

            with self.assertRaisesRegex(ValueError, "empty or non-finite"):
                load_percolator_psms(target_path, decoy_path)

    def test_global_tsv_can_be_filtered_by_run_before_duplicate_check(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            target_path = root / "global_target.tsv"
            decoy_path = root / "global_decoy.tsv"
            pd.DataFrame(
                [
                    ["run_a::run.101.101.2_1", 5.2, 0.001, 0.002, "K.PEPTIDE.R", "P1"],
                    ["run_b::run.101.101.2_1", 4.8, 0.002, 0.003, "K.PEPTIDE.R", "P1"],
                ],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(target_path, sep="\t", index=False)
            pd.DataFrame(columns=PERCOLATOR_COLUMNS).to_csv(
                decoy_path, sep="\t", index=False
            )

            psms = load_percolator_psms(
                target_path, decoy_path, run_id="run_a"
            )

            self.assertEqual(psms["psm_id"].tolist(), ["101_2_PEPTIDE"])
            self.assertEqual(psms["score"].tolist(), [5.2])

    def test_fdr_threshold_filters_targets_and_keeps_decoys(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_path = root / "input.parquet"
            output_path = root / "output.parquet"
            target_path = root / "target.tsv"
            decoy_path = root / "decoy.tsv"
            pd.DataFrame(
                {
                    "scan": [1, 2, 3],
                    "charge": [2, 2, 2],
                    "precursor_sequence": ["GOODSEQ", "BADSEQ", "DECOYSEQ"],
                    "label": [1, 1, 0],
                }
            ).to_parquet(input_path, index=False)
            pd.DataFrame(
                [
                    ["run.1.1.2_1", 5.0, 0.01, 0.01, "K.GOODSEQ.R", "P1"],
                    ["run.2.2.2_1", 4.0, 0.20, 0.20, "K.BADSEQ.R", "P2"],
                ],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(target_path, sep="\t", index=False)
            pd.DataFrame(
                [["run.3.3.2_1", -1.0, 0.90, 0.90, "K.DECOYSEQ.R", "rev_P3"]],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(decoy_path, sep="\t", index=False)

            enrich_parquet_with_percolator(
                input_path,
                target_path,
                decoy_path,
                output_path,
                fdr_threshold=0.05,
            )

            result = pd.read_parquet(output_path)
            self.assertEqual(result["scan"].tolist(), [1, 3])


if __name__ == "__main__":
    unittest.main()
