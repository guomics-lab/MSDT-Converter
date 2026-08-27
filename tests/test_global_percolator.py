import subprocess
import tempfile
import unittest
from pathlib import Path

from scripts.percolator import run_global_percolator


PIN_HEADER = (
    "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\thyperscore\tPeptide\tProteins\n"
)


class GlobalPercolatorTests(unittest.TestCase):
    def test_compatible_pins_are_combined_and_scored_once(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            pin_a = root / "a.pin"
            pin_b = root / "b.pin"
            pin_a.write_text(
                PIN_HEADER
                + "DefaultDirection\t-\t-\t-\t-\t1\t-\t-\n"
                + "run.101.101.2_1\t1\t101\t1000\t999\t40\tK.PEPTIDE.R\tP1\n",
                encoding="utf-8",
            )
            pin_b.write_text(
                PIN_HEADER
                + "DefaultDirection\t-\t-\t-\t-\t1\t-\t-\n"
                + "run.101.101.2_1\t-1\t101\t1000\t999\t10\tK.PEPTIDE.R\trev_P1\n",
                encoding="utf-8",
            )
            percolator_exe = root / "percolator"
            percolator_exe.write_text("test executable", encoding="utf-8")
            calls = []

            def fake_runner(command, **kwargs):
                calls.append((command, kwargs))
                target = Path(command[command.index("--results-psms") + 1])
                decoy = Path(command[command.index("--decoy-results-psms") + 1])
                target.write_text("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\n", encoding="utf-8")
                decoy.write_text("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, "", "")

            artifacts = run_global_percolator(
                {"run_a": pin_a, "run_b": pin_b},
                percolator_exe,
                root / "global",
                threads=4,
                runner=fake_runner,
            )

            combined_lines = artifacts.combined_pin.read_text(
                encoding="utf-8"
            ).splitlines()
            self.assertEqual(combined_lines.count(PIN_HEADER.rstrip("\n")), 1)
            self.assertEqual(
                sum(line.startswith("DefaultDirection\t") for line in combined_lines),
                1,
            )
            self.assertIn("run_a::run.101.101.2_1", combined_lines[2])
            self.assertIn("run_b::run.101.101.2_1", combined_lines[3])
            header = combined_lines[0].split("\t")
            scan_index = header.index("ScanNr")
            run_a_scan = combined_lines[2].split("\t")[scan_index]
            run_b_scan = combined_lines[3].split("\t")[scan_index]
            self.assertNotEqual(run_a_scan, run_b_scan)
            self.assertEqual(len(calls), 1)
            self.assertIn("--post-processing-tdc", calls[0][0])
            self.assertIn("4", calls[0][0])
            self.assertTrue(artifacts.target_tsv.is_file())
            self.assertTrue(artifacts.decoy_tsv.is_file())

    def test_candidates_from_one_spectrum_share_the_remapped_scan_number(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            pin = root / "a.pin"
            pin.write_text(
                PIN_HEADER
                + "run.101.101.2_1\t1\t101\t1000\t999\t40\tK.PEPTIDE.R\tP1\n"
                + "run.101.101.3_1\t1\t101\t1001\t998\t30\tK.OTHERSEQ.R\tP2\n",
                encoding="utf-8",
            )
            percolator_exe = root / "percolator"
            percolator_exe.touch()

            def fake_runner(command, **kwargs):
                Path(command[command.index("--results-psms") + 1]).touch()
                Path(command[command.index("--decoy-results-psms") + 1]).touch()

            artifacts = run_global_percolator(
                {"run_a": pin}, percolator_exe, root / "global", runner=fake_runner
            )

            lines = artifacts.combined_pin.read_text(encoding="utf-8").splitlines()
            scan_index = lines[0].split("\t").index("ScanNr")
            self.assertEqual(
                lines[1].split("\t")[scan_index],
                lines[2].split("\t")[scan_index],
            )


if __name__ == "__main__":
    unittest.main()
