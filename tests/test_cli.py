import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from convert import _resolve_global_run_id, main


PERCOLATOR_COLUMNS = [
    "PSMId",
    "score",
    "q-value",
    "posterior_error_prob",
    "peptide",
    "proteinIds",
]


class CliTests(unittest.TestCase):
    def test_fp_msdt_accepts_pin_without_optional_msbooster_columns(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            raw = root / "raw.parquet"
            pin = root / "sample.pin"
            target = root / "target.tsv"
            decoy = root / "decoy.tsv"
            output = root / "output.parquet"
            pd.DataFrame(
                {
                    "scan": [101],
                    "precursor_mz": [500.2],
                    "rt": [10.0],
                    "mz_array": [[100.0, 200.0]],
                    "intensity_array": [[10.0, 20.0]],
                }
            ).to_parquet(raw, index=False)
            pd.DataFrame(
                [["run.101.101.2_1", 1, 101, 1000.0, 10.0, 1, 0, 40.0,
                  5.0, 8, "b,y", "K.PEPTIDE.R", "P1"]],
                columns=[
                    "SpecId", "Label", "ScanNr", "ExpMass", "retentiontime",
                    "rank", "isotope_errors", "hyperscore", "delta_hyperscore",
                    "matched_ion_num", "ion_series", "Peptide", "Proteins",
                ],
            ).to_csv(pin, sep="\t", index=False)
            pd.DataFrame(
                [["run.101.101.2_1", 5.0, 0.001, 0.002, "K.PEPTIDE.R", "P1"]],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(target, sep="\t", index=False)
            pd.DataFrame(columns=PERCOLATOR_COLUMNS).to_csv(
                decoy, sep="\t", index=False
            )

            exit_code = main(
                [
                    "fp-msdt", "--instrument", "mzml",
                    "--raw-spectrum", str(raw), "--pin", str(pin),
                    "--target-tsv", str(target), "--decoy-tsv", str(decoy),
                    "--output", str(output),
                ]
            )

            result = pd.read_parquet(output)
            self.assertEqual(exit_code, 0)
            self.assertTrue(result["unweighted_spectral_entropy"].isna().all())
            self.assertTrue(result["delta_RT_loess"].isna().all())
            self.assertEqual(result["score"].tolist(), [5.0])

    def test_enrich_command_adds_percolator_fields(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            parquet = root / "input.parquet"
            target = root / "target.tsv"
            decoy = root / "decoy.tsv"
            output = root / "output.parquet"
            pd.DataFrame(
                {
                    "scan": [101],
                    "charge": [2],
                    "precursor_sequence": ["PEPTIDE"],
                    "label": [1],
                }
            ).to_parquet(parquet, index=False)
            pd.DataFrame(
                [["run.101.101.2_1", 5.0, 0.01, 0.02, "K.PEPTIDE.R", "P1"]],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(target, sep="\t", index=False)
            pd.DataFrame(columns=PERCOLATOR_COLUMNS).to_csv(
                decoy, sep="\t", index=False
            )

            exit_code = main(
                [
                    "enrich",
                    "--parquet",
                    str(parquet),
                    "--target-tsv",
                    str(target),
                    "--decoy-tsv",
                    str(decoy),
                    "--output",
                    str(output),
                    "--global-fdr",
                    "0.01",
                ]
            )

            result = pd.read_parquet(output)
            self.assertEqual(exit_code, 0)
            self.assertEqual(result["score"].tolist(), [5.0])
            self.assertIn("q-value", result.columns)
            self.assertIn("PEP", result.columns)

    def test_legacy_config_invocation_remains_valid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            config = Path(temp_dir) / "empty.json"
            config.write_text(json.dumps({}), encoding="utf-8")

            self.assertEqual(main(["-config", str(config)]), 0)

    def test_global_run_id_comes_from_pin_mapping_key(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            pin = Path(temp_dir) / "filename_does_not_match_run.pin"
            pin.write_text("fixture", encoding="utf-8")

            self.assertEqual(
                _resolve_global_run_id(pin, {"declared_run": str(pin)}),
                "declared_run",
            )


if __name__ == "__main__":
    unittest.main()
