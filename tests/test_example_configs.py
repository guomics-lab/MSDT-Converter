import json
import unittest
from pathlib import Path, PurePosixPath


REPO_ROOT = Path(__file__).resolve().parents[1]


class ExampleConfigTests(unittest.TestCase):
    def test_packaged_default_workflow_keeps_percolator_tsv_files(self):
        workflow = REPO_ROOT / "workflows" / "Default-v2.workflow"

        self.assertTrue(workflow.is_file())
        settings = {
            line.split("=", 1)[0]: line.split("=", 1)[1]
            for line in workflow.read_text(encoding="utf-8").splitlines()
            if "=" in line
        }
        self.assertEqual(settings["percolator.run-percolator"], "true")
        self.assertEqual(settings["percolator.keep-tsv-files"], "true")

    def test_per_run_fragpipe_outputs_share_one_directory_and_have_no_run_id(self):
        examples = (
            ("config_mzml.json", "mzml"),
            ("config_wiff.json", "wiff"),
        )

        for config_name, data_type in examples:
            with self.subTest(config=config_name, data_type=data_type):
                config = json.loads(
                    (REPO_ROOT / config_name).read_text(encoding="utf-8")
                )
                self.assertEqual(
                    config["generate_fragpipe_search_result"]["workflow_path"],
                    "workflows/Default-v2.workflow",
                )
                params = config["generate_msdt"][data_type]
                pin_parent = PurePosixPath(params["fp_pin_path"]).parent
                target_parent = PurePosixPath(
                    params["percolator_target_path"]
                ).parent
                decoy_parent = PurePosixPath(
                    params["percolator_decoy_path"]
                ).parent

                self.assertEqual(target_parent, pin_parent)
                self.assertEqual(decoy_parent, pin_parent)
                self.assertNotIn("run_id", params)


if __name__ == "__main__":
    unittest.main()
