import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.search_engine import (
    build_manifest,
    generate_fp_search_result_fn,
    prepare_workflow,
    read_file_list,
    run_fragpipe,
)


def write_complete_fragpipe_outputs(result_dir, stem, *, edited=True):
    pin_name = f"{stem}_edited.pin" if edited else f"{stem}.pin"
    (result_dir / pin_name).write_text(
        "SpecId\tLabel\n", encoding="utf-8"
    )
    for kind in ("target", "decoy"):
        (result_dir / f"{stem}_percolator_{kind}_psms.tsv").write_text(
            "PSMId\tscore\n", encoding="utf-8"
        )


class FragPipeBatchTests(unittest.TestCase):
    def test_fragpipe_uses_decoy_fasta_created_by_current_run(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            fragpipe = root / "FragPipe-21.1" / "bin" / "fragpipe"
            fragpipe.parent.mkdir(parents=True)
            fragpipe.touch()
            fasta = root / "database" / "database.fasta"
            fasta.parent.mkdir()
            fasta.write_text(">P1\nPEPTIDE\n", encoding="utf-8")
            stale_decoy = fasta.parent / "aaa_old_database.fasta.fas"
            stale_decoy.write_text(">rev_old\nOLD\n", encoding="utf-8")
            fresh_decoy = fasta.parent / "zzz_new_database.fasta.fas"
            workflow = root / "default.workflow"
            workflow.write_text("database.db-path=/old/path\n", encoding="utf-8")
            manifest = root / "input.fp-manifest"
            manifest.touch()
            output = root / "results"
            output.mkdir()

            def fake_philosopher(command, **kwargs):
                if "database" in command and "--custom" in command:
                    fresh_decoy.write_text(">rev_new\nNEW\n", encoding="utf-8")

            with patch("scripts.search_engine.subprocess.run", fake_philosopher), patch(
                "scripts.search_engine.run_cmd"
            ):
                run_fragpipe(
                    manifest,
                    workflow,
                    output,
                    fragpipe,
                    2,
                    fasta,
                )

            generated_workflow = output / workflow.name
            self.assertIn(
                f"database.db-path={fresh_decoy}",
                generated_workflow.read_text(encoding="utf-8"),
            )
            self.assertNotIn(
                f"database.db-path={stale_decoy}",
                generated_workflow.read_text(encoding="utf-8"),
            )

    def test_one_column_file_list_becomes_batch_manifest(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            sample_a = root / "sample a.mzML"
            sample_b = root / "sample_b.mzML"
            sample_a.touch()
            sample_b.touch()
            file_list = root / "file_list.txt"
            file_list.write_text(
                f"# inputs\n{sample_a}\n{sample_b}\n", encoding="utf-8"
            )

            entries = read_file_list(file_list)
            manifest = build_manifest(entries, root / "results")

            self.assertEqual(len(entries), 2)
            self.assertEqual(
                manifest.read_text(encoding="utf-8").splitlines(),
                [
                    f"{sample_a}\texp\t\tDDA",
                    f"{sample_b}\texp\t\tDDA",
                ],
            )

    def test_official_four_column_manifest_is_preserved(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            sample = root / "sample.d"
            sample.mkdir()
            file_list = root / "file_list.tsv"
            file_list.write_text(
                f"{sample}\tcontrol\t2\tDDA\n", encoding="utf-8"
            )

            entries = read_file_list(file_list)

            self.assertEqual(entries[0].experiment, "control")
            self.assertEqual(entries[0].bioreplicate, "2")
            self.assertEqual(entries[0].data_type, "DDA")

    def test_workflow_copy_keeps_percolator_tsv_and_updates_database(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "default.workflow"
            destination = root / "run.workflow"
            source.write_text(
                "database.db-path=/old/db.fasta\n"
                "percolator.keep-tsv-files=false\n"
                "percolator.run-percolator=true\n",
                encoding="utf-8",
            )

            prepare_workflow(source, destination, root / "decoy.fasta")

            self.assertIn(
                f"database.db-path={root / 'decoy.fasta'}",
                destination.read_text(encoding="utf-8"),
            )
            self.assertIn(
                "percolator.keep-tsv-files=true",
                destination.read_text(encoding="utf-8"),
            )
            self.assertIn(
                "percolator.keep-tsv-files=false",
                source.read_text(encoding="utf-8"),
            )

    def test_fragpipe_search_reads_all_inputs_from_file_list(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            samples = [root / "sample_a.mzML", root / "sample_b.mzML"]
            for sample in samples:
                sample.touch()
            file_list = root / "file_list.txt"
            file_list.write_text(
                "\n".join(str(sample) for sample in samples) + "\n",
                encoding="utf-8",
            )
            workflow = root / "default.workflow"
            workflow.write_text("percolator.run-percolator=true\n", encoding="utf-8")
            fasta = root / "database.fasta"
            fasta.write_text(">P1\nPEPTIDE\n", encoding="utf-8")
            workdir = root / "results"
            observed_manifest_lines = []

            def fake_run_fragpipe(manifest_path, *args):
                observed_manifest_lines.extend(
                    Path(manifest_path).read_text(encoding="utf-8").splitlines()
                )
                result_dir = Path(args[1]) / "exp"
                result_dir.mkdir(parents=True)
                for sample in samples:
                    write_complete_fragpipe_outputs(result_dir, sample.stem)

            with patch("scripts.search_engine.run_fragpipe", fake_run_fragpipe):
                state = generate_fp_search_result_fn(
                    {
                        "file_list": str(file_list),
                        "workdir": str(workdir),
                        "fasta_path": str(fasta),
                        "workflow_path": str(workflow),
                        "thread_num": 4,
                    }
                )

            self.assertEqual(state, 0)
            self.assertEqual(len(observed_manifest_lines), 2)

    def test_plain_pin_with_percolator_tsvs_is_a_complete_search(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            sample = root / "sample.mzML"
            sample.touch()
            file_list = root / "file_list.txt"
            file_list.write_text(f"{sample}\n", encoding="utf-8")
            workflow = root / "default.workflow"
            workflow.write_text("percolator.run-percolator=true\n", encoding="utf-8")
            fasta = root / "database.fasta"
            fasta.write_text(">P1\nPEPTIDE\n", encoding="utf-8")
            workdir = root / "results"
            result_dir = workdir / "exp"
            result_dir.mkdir(parents=True)
            write_complete_fragpipe_outputs(result_dir, "sample", edited=False)

            with patch("scripts.search_engine.run_fragpipe") as run_fragpipe_mock:
                state = generate_fp_search_result_fn(
                    {
                        "file_list": str(file_list),
                        "workdir": str(workdir),
                        "fasta_path": str(fasta),
                        "workflow_path": str(workflow),
                        "thread_num": 2,
                    }
                )

            self.assertEqual(state, 1)
            run_fragpipe_mock.assert_not_called()

    def test_zero_byte_percolator_outputs_are_not_reported_complete(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            sample = root / "sample.mzML"
            sample.touch()
            file_list = root / "file_list.txt"
            file_list.write_text(f"{sample}\n", encoding="utf-8")
            workflow = root / "default.workflow"
            workflow.touch()
            fasta = root / "database.fasta"
            fasta.touch()
            result_dir = root / "results" / "exp"
            result_dir.mkdir(parents=True)
            (result_dir / "sample_edited.pin").write_text(
                "SpecId\tLabel\n", encoding="utf-8"
            )
            (result_dir / "sample_percolator_target_psms.tsv").touch()
            (result_dir / "sample_percolator_decoy_psms.tsv").touch()

            with patch("scripts.search_engine.run_fragpipe") as run_fragpipe_mock:
                state = generate_fp_search_result_fn(
                    {
                        "file_list": str(file_list),
                        "workdir": str(root / "results"),
                        "fasta_path": str(fasta),
                        "workflow_path": str(workflow),
                        "thread_num": 2,
                    }
                )

            self.assertEqual(state, -1)
            run_fragpipe_mock.assert_called_once()

    def test_pin_without_percolator_tsvs_is_not_reported_complete(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            sample = root / "sample.mzML"
            sample.touch()
            file_list = root / "file_list.txt"
            file_list.write_text(f"{sample}\n", encoding="utf-8")
            workflow = root / "default.workflow"
            workflow.touch()
            fasta = root / "database.fasta"
            fasta.touch()
            workdir = root / "results"
            workdir.mkdir()
            (workdir / "sample_edited.pin").touch()

            with patch("scripts.search_engine.run_fragpipe"):
                state = generate_fp_search_result_fn(
                    {
                        "file_list": str(file_list),
                        "workdir": str(workdir),
                        "fasta_path": str(fasta),
                        "workflow_path": str(workflow),
                        "thread_num": 2,
                    }
                )

            self.assertEqual(state, -1)

    def test_outputs_from_different_directories_do_not_form_a_result_set(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            sample = root / "sample.mzML"
            sample.touch()
            file_list = root / "file_list.tsv"
            file_list.write_text(f"{sample}\texp\t\tDDA\n", encoding="utf-8")
            workflow = root / "default.workflow"
            workflow.touch()
            fasta = root / "database.fasta"
            fasta.touch()
            workdir = root / "results"
            pin_dir = workdir / "exp" / "old_pin"
            tsv_dir = workdir / "exp" / "new_tsv"
            pin_dir.mkdir(parents=True)
            tsv_dir.mkdir(parents=True)
            (pin_dir / "sample.pin").touch()
            (tsv_dir / "sample_percolator_target_psms.tsv").touch()
            (tsv_dir / "sample_percolator_decoy_psms.tsv").touch()

            with patch("scripts.search_engine.run_fragpipe"):
                state = generate_fp_search_result_fn(
                    {
                        "file_list": str(file_list), "workdir": str(workdir),
                        "fasta_path": str(fasta), "workflow_path": str(workflow),
                        "thread_num": 2,
                    }
                )

            self.assertEqual(state, -1)

    def test_same_input_stem_is_resolved_by_manifest_experiment(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            samples = [root / "a" / "sample.mzML", root / "b" / "sample.mzML"]
            for sample in samples:
                sample.parent.mkdir()
                sample.touch()
            file_list = root / "file_list.tsv"
            file_list.write_text(
                f"{samples[0]}\texp_a\t\tDDA\n"
                f"{samples[1]}\texp_b\t\tDDA\n",
                encoding="utf-8",
            )
            workflow = root / "default.workflow"
            workflow.touch()
            fasta = root / "database.fasta"
            fasta.touch()
            workdir = root / "results"
            for experiment in ("exp_a", "exp_b"):
                result_dir = workdir / experiment
                result_dir.mkdir(parents=True)
                write_complete_fragpipe_outputs(result_dir, "sample")

            with patch("scripts.search_engine.run_fragpipe") as run_fragpipe_mock:
                state = generate_fp_search_result_fn(
                    {
                        "file_list": str(file_list), "workdir": str(workdir),
                        "fasta_path": str(fasta), "workflow_path": str(workflow),
                        "thread_num": 2,
                    }
                )

            self.assertEqual(state, 1)
            run_fragpipe_mock.assert_not_called()


if __name__ == "__main__":
    unittest.main()
