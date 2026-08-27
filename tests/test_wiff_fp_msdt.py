import subprocess
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.generate_msdt import gen_wiff_fragpipe_msdt


PIN_COLUMNS = [
    "SpecId",
    "Label",
    "ScanNr",
    "ExpMass",
    "retentiontime",
    "rank",
    "isotope_errors",
    "hyperscore",
    "delta_hyperscore",
    "matched_ion_num",
    "ion_series",
    "unweighted_spectral_entropy",
    "delta_RT_loess",
    "Peptide",
    "Proteins",
]

PERCOLATOR_COLUMNS = [
    "PSMId",
    "score",
    "q-value",
    "posterior_error_prob",
    "peptide",
    "proteinIds",
]


class WiffFragPipeMsdtTests(unittest.TestCase):
    def test_search_scans_join_to_native_scans_and_multiple_psms_are_kept(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            native_raw = root / "native_raw.parquet"
            wiff_mzml = root / "sample.mzML"
            fp_pin = root / "sample_edited.pin"
            target_tsv = root / "sample_percolator_target_psms.tsv"
            decoy_tsv = root / "sample_percolator_decoy_psms.tsv"
            output = root / "sample_fp_wiff_msdt.parquet"
            wiff_mzml.write_text("fixture", encoding="utf-8")
            pd.DataFrame(
                {
                    "scan": [100, 101],
                    "precursor_mz": [500.2, 600.3],
                    "rt": [10.0, 11.0],
                    "mz_array": ["100.0,200.0", "110.0,210.0"],
                    "intensity_array": ["10.0,20.0", "11.0,21.0"],
                }
            ).to_parquet(native_raw, index=False)
            pd.DataFrame(
                [
                    ["run.101.101.2_1", 1, 101, 1000, 10.0, 1, 0, 40, 5, 8, "b,y", 0.8, 0.1, "K.PEPTIDE.R", "P1"],
                    ["run.101.101.3_1", 1, 101, 1001, 10.0, 1, 0, 35, 4, 7, "b,y", 0.7, 0.2, "K.OTHERSEQ.R", "P2"],
                    ["run.102.102.2_1", -1, 102, 1200, 11.0, 1, 0, 10, 1, 3, "b,y", 0.2, 0.3, "K.DECOYSEQ.R", "rev_P3"],
                ],
                columns=PIN_COLUMNS,
            ).to_csv(fp_pin, sep="\t", index=False)
            pd.DataFrame(
                [
                    ["run.101.101.2_1", 5.0, 0.001, 0.002, "K.PEPTIDE.R", "P1"],
                    ["run.101.101.3_1", 4.0, 0.002, 0.003, "K.OTHERSEQ.R", "P2"],
                ],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(target_tsv, sep="\t", index=False)
            pd.DataFrame(
                [["run.102.102.2_1", -1.0, 0.5, 0.8, "K.DECOYSEQ.R", "rev_P3"]],
                columns=PERCOLATOR_COLUMNS,
            ).to_csv(decoy_tsv, sep="\t", index=False)

            def fake_runner(command, **kwargs):
                search_scan_tsv = Path(command[2])
                pd.DataFrame(
                    {
                        "scan": [
                            "sample=1 period=1 cycle=101 experiment=2",
                            "sample=1 period=1 cycle=101 experiment=3",
                        ],
                        "precursor_mz": [500.2, 600.3],
                        "rt": [10.0, 11.0],
                    }
                ).to_csv(
                    search_scan_tsv, sep="\t", index=False
                )
                return subprocess.CompletedProcess(command, 0, "", "")

            state = gen_wiff_fragpipe_msdt(
                native_raw,
                wiff_mzml,
                fp_pin,
                target_tsv,
                decoy_tsv,
                output,
                unify_residue=True,
                mzml_extractor=root / "linux_mzml_rawspectrum",
                runner=fake_runner,
            )

            result = pd.read_parquet(output)
            self.assertEqual(state, 0)
            self.assertEqual(result["scan"].tolist(), [100, 100, 101])
            self.assertEqual(result["charge"].tolist(), [2, 3, 2])
            self.assertEqual(result["score"].tolist(), [5.0, 4.0, -1.0])
            self.assertNotIn("search_scan", result.columns)
            self.assertFalse(result[["score", "q-value", "PEP"]].isna().any().any())

    def test_reordered_native_and_search_spectra_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            native_raw = root / "native_raw.parquet"
            wiff_mzml = root / "sample.mzML"
            fp_pin = root / "sample_edited.pin"
            target_tsv = root / "target.tsv"
            decoy_tsv = root / "decoy.tsv"
            output = root / "output.parquet"
            wiff_mzml.write_text("fixture", encoding="utf-8")
            pd.DataFrame(
                {
                    "scan": [1, 2],
                    "precursor_mz": [500.2, 600.3],
                    "rt": [10.0, 11.0],
                    "mz_array": ["100.0", "110.0"],
                    "intensity_array": ["10.0", "11.0"],
                }
            ).to_parquet(native_raw, index=False)
            pd.DataFrame(columns=PIN_COLUMNS).to_csv(fp_pin, sep="\t", index=False)
            pd.DataFrame(columns=PERCOLATOR_COLUMNS).to_csv(
                target_tsv, sep="\t", index=False
            )
            pd.DataFrame(columns=PERCOLATOR_COLUMNS).to_csv(
                decoy_tsv, sep="\t", index=False
            )

            def fake_runner(command, **kwargs):
                pd.DataFrame(
                    {
                        "scan": [102, 101],
                        "precursor_mz": [600.3, 500.2],
                        "rt": [11.0, 10.0],
                    }
                ).to_csv(Path(command[2]), sep="\t", index=False)
                return subprocess.CompletedProcess(command, 0, "", "")

            state = gen_wiff_fragpipe_msdt(
                native_raw,
                wiff_mzml,
                fp_pin,
                target_tsv,
                decoy_tsv,
                output,
                mzml_extractor=root / "linux_mzml_rawspectrum",
                runner=fake_runner,
            )

            self.assertEqual(state, -1)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
