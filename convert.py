"""MSDT-Converter command-line entry point."""

from __future__ import annotations

import argparse
import copy
import json
import logging
import os.path
from pathlib import Path
from typing import Sequence

from scripts.generate_rawspectrum import generate_rawspectrum_fn
from scripts.generate_msdt import (
    gen_mzml_fragpipe_msdt,
    gen_wiff_fragpipe_msdt,
    generate_msdt_fn,
)
from scripts.mgf2parquet import mgf_to_parquet
from scripts.percolator import (
    enrich_parquet_with_percolator,
    run_global_percolator,
)
from scripts.search_engine import (
    generate_fp_search_result_fn,
    generate_sage_search_result_fn,
)
from scripts.msdt2mgf import msdt2mgf


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def load_config(config_path: str | Path) -> dict:
    """Read JSON configuration and recursively remove comment fields."""
    source = Path(config_path)
    if not source.is_file():
        raise FileNotFoundError(f"Configuration file not found: {source}")
    with source.open("r", encoding="utf-8") as handle:
        config = json.load(handle)

    def remove_comments(value):
        if isinstance(value, dict):
            return {
                key: remove_comments(item)
                for key, item in value.items()
                if not key.startswith("_comment")
            }
        if isinstance(value, list):
            return [remove_comments(item) for item in value]
        return value

    return remove_comments(config)


def parse_config(config: dict) -> dict:
    """Select enabled pipeline steps while preserving their parameters."""
    steps = {}
    raw = config.get("generate_rawspectrum", {})
    if raw.get("need", False):
        steps["generate_rawspectrum"] = {
            "data_type": raw["data_type"],
            "input": raw["data_path"],
            "output": raw["output"],
        }
    sage = config.get("generate_sage_search_result", {})
    if sage.get("need", False):
        steps["generate_sage_search_result"] = sage
    fragpipe = config.get("generate_fragpipe_search_result", {})
    if fragpipe.get("need", False):
        steps["generate_fragpipe_search_result"] = fragpipe
    global_percolator = config.get("global_percolator", {})
    if global_percolator.get("need", False):
        steps["global_percolator"] = global_percolator
    msdt = config.get("generate_msdt", {})
    if msdt.get("need", False):
        steps["generate_msdt"] = {
            "tims": msdt.get("tims", {}),
            "mzml": msdt.get("mzml", {}),
            "wiff": msdt.get("wiff", {}),
        }
    mgf = config.get("convert_2_msdt", {}).get("mgf", {})
    if mgf.get("need", False):
        steps["convert_2_msdt"] = {"mgf": mgf}
    to_mgf = config.get("msdt_2_mgf", {})
    if to_mgf.get("need", False):
        steps["msdt2mgf"] = to_mgf
    return steps


def _resolve_global_run_id(
    pin_path: str | Path, pin_files: dict[str, str | Path]
) -> str:
    """Resolve a PIN to the exact run ID used to prefix the global PIN."""
    requested = Path(pin_path).expanduser().resolve(strict=False)
    matches = [
        run_id
        for run_id, candidate in pin_files.items()
        if Path(candidate).expanduser().resolve(strict=False) == requested
    ]
    if len(matches) != 1:
        raise ValueError(
            f"PIN path {pin_path} must occur exactly once in "
            "global_percolator.pin_files"
        )
    return matches[0]


def execute_steps(steps: dict) -> int:
    """Execute configured steps and return a process-style exit code."""
    states = []
    if not steps:
        logger.info("No steps to execute.")
        return 0
    if "generate_rawspectrum" in steps:
        states.append(generate_rawspectrum_fn(steps["generate_rawspectrum"]))
    if "generate_sage_search_result" in steps:
        states.append(
            generate_sage_search_result_fn(steps["generate_sage_search_result"])
        )
    if "generate_fragpipe_search_result" in steps:
        states.append(
            generate_fp_search_result_fn(steps["generate_fragpipe_search_result"])
        )

    global_artifacts = None
    global_options = steps.get("global_percolator")
    if global_options:
        pin_files = global_options.get("pin_files")
        if not isinstance(pin_files, dict) or not pin_files:
            raise ValueError(
                "global_percolator.pin_files must map run_id to PIN path"
            )
        global_artifacts = run_global_percolator(
            pin_files,
            global_options["percolator_executable"],
            global_options["output_dir"],
            threads=int(global_options.get("threads", 1)),
        )

    if "generate_msdt" in steps:
        msdt_params = copy.deepcopy(steps["generate_msdt"])
        if global_artifacts is not None:
            threshold = global_options.get("fdr_threshold")
            for data_type in ("mzml", "wiff"):
                params = msdt_params.get(data_type, {})
                if params.get("need_fragpipe", False):
                    params["percolator_target_path"] = str(
                        global_artifacts.target_tsv
                    )
                    params["percolator_decoy_path"] = str(
                        global_artifacts.decoy_tsv
                    )
                    params["fdr_threshold"] = threshold
                    if not params.get("fp_pin_path"):
                        raise ValueError(
                            f"generate_msdt.{data_type}.fp_pin_path is required "
                            "for global FDR"
                        )
                    mapped_run_id = _resolve_global_run_id(
                        params["fp_pin_path"], pin_files
                    )
                    configured_run_id = params.get("run_id")
                    if configured_run_id and configured_run_id != mapped_run_id:
                        raise ValueError(
                            f"generate_msdt.{data_type}.run_id "
                            f"({configured_run_id}) does not match the global PIN "
                            f"mapping key ({mapped_run_id})"
                        )
                    params["run_id"] = mapped_run_id
        states.append(generate_msdt_fn(msdt_params))
    if "convert_2_msdt" in steps:
        states.append(mgf_to_parquet(steps["convert_2_msdt"]["mgf"]))
    if "msdt2mgf" in steps:
        states.append(msdt2mgf(steps["msdt2mgf"]))
    return 0 if all(state in (0, 1) for state in states) else 1


def run_config(
    config_path: str | Path,
    *,
    file_list: str | None = None,
    threads: int | None = None,
    global_fdr: float | None = None,
    percolator_executable: str | None = None,
) -> int:
    """Run a JSON workflow with optional command-line overrides."""
    config = load_config(config_path)
    fragpipe = config.get("generate_fragpipe_search_result", {})
    if file_list is not None:
        fragpipe["file_list"] = file_list
        fragpipe.pop("data_path", None)
    if threads is not None:
        fragpipe["thread_num"] = threads
    if global_fdr is not None:
        global_options = config.get("global_percolator")
        if not global_options:
            raise ValueError(
                "--global-fdr requires a global_percolator section in the config"
            )
        global_options["need"] = True
        global_options["fdr_threshold"] = global_fdr
        if percolator_executable is not None:
            global_options["percolator_executable"] = percolator_executable
    return execute_steps(parse_config(config))


def _pin_mapping(values: list[str]) -> dict[str, str]:
    mapping = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--pin must use RUN_ID=PATH format")
        run_id, path = value.split("=", 1)
        if not run_id or not path or run_id in mapping:
            raise ValueError(f"Invalid or duplicate --pin value: {value}")
        mapping[run_id] = path
    return mapping


def _add_global_fdr_argument(command_parser: argparse.ArgumentParser) -> None:
    command_parser.add_argument(
        "--global-fdr",
        type=float,
        help="target q-value threshold; all decoy PSMs are retained",
    )


def _add_run_id_argument(command_parser: argparse.ArgumentParser) -> None:
    command_parser.add_argument(
        "--run-id",
        help=(
            "run key used only with pooled global-Percolator TSVs whose "
            "PSMId values start with RUN_ID::; omit for per-run TSVs"
        ),
    )


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="MSDT-Converter v2")
    parser.add_argument(
        "-config", "--config", dest="legacy_config", help="legacy config JSON"
    )
    commands = parser.add_subparsers(dest="command")

    run = commands.add_parser("run", help="run a JSON pipeline")
    run.add_argument("--config", required=True)
    run.add_argument("--file-list")
    run.add_argument("--threads", type=int)
    _add_global_fdr_argument(run)
    run.add_argument("--percolator-exe")

    search = commands.add_parser("fp-search", help="run batch FragPipe search")
    search.add_argument("--file-list", required=True)
    search.add_argument("--workdir", required=True)
    search.add_argument("--fasta", required=True)
    search.add_argument("--workflow", required=True)
    search.add_argument("--threads", type=int, default=1)

    enrich = commands.add_parser(
        "enrich", help="add Percolator fields to an FP MSDT Parquet"
    )
    enrich.add_argument("--file", required=True)
    enrich.add_argument("--data_type", required=True)
    enrich.add_argument("--workdir", required=True)
    enrich.add_argument("--fasta", required=True)
    enrich.add_argument("--workflow", required=True)
    enrich.add_argument("--threads", type=int, default=1)
    _add_run_id_argument(enrich)
    _add_global_fdr_argument(enrich)

    build = commands.add_parser("fp-msdt", help="build an FP-derived MSDT")
    build.add_argument("--instrument", choices=("mzml", "wiff"), required=True)
    build.add_argument("--raw-spectrum", required=True)
    build.add_argument("--pin", required=True)
    build.add_argument("--target-tsv", required=True)
    build.add_argument("--decoy-tsv", required=True)
    build.add_argument("--output", required=True)
    build.add_argument(
        "--wiff-mzml",
        help=(
            "mzML converted from the SCIEX .wiff/.wiff.scan pair; "
            "native WIFF conversion is performed outside this Linux tool"
        ),
    )
    _add_run_id_argument(build)
    _add_global_fdr_argument(build)
    build.add_argument("--no-unify-residue", action="store_true")

    global_command = commands.add_parser(
        "global-percolator",
        help="pool PINs, remap cross-run ScanNr values, and run Percolator once",
    )
    global_command.add_argument(
        "--pin", action="append", required=True, metavar="RUN_ID=PATH"
    )
    global_command.add_argument("--percolator-exe", required=True)
    global_command.add_argument("--output-dir", required=True)
    global_command.add_argument("--threads", type=int, default=1)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            return run_config(
                args.config,
                file_list=args.file_list,
                threads=args.threads,
                global_fdr=args.global_fdr,
                percolator_executable=args.percolator_exe,
            )
        if args.command == "fp-search":
            state = generate_fp_search_result_fn(
                {
                    "file_list": args.file_list,
                    "workdir": args.workdir,
                    "fasta_path": args.fasta,
                    "workflow_path": args.workflow,
                    "thread_num": args.threads,
                }
            )
            return 0 if state in (0, 1) else 1
        if args.command == "enrich":

            base_file_name = os.path.basename(args.file).removesuffix('.mzML')
            rawspectrum_parquet = os.path.join(args.workdir, base_file_name + "_rawspectrum.parquet")
            target_tsv = os.path.join(args.workdir, 'exp', base_file_name + "_percolator_target_psms.tsv")
            decoy_tsv = os.path.join(args.workdir, 'exp', base_file_name + "_percolator_decoy_psms.tsv")
            final_output_file = os.path.join(args.workdir, base_file_name + "_fp_msdt_v2.parquet")
            pin_file = os.path.join(args.workdir, 'exp', base_file_name + "_edited.pin")

            generate_rawspectrum_fn_param = {
                "data_type": args.data_type,
                "input": args.file,
                "output": rawspectrum_parquet,
            }
            generate_rawspectrum_fn(generate_rawspectrum_fn_param)

            generate_fp_search_result_fn(
                {
                    "data_path": args.file,
                    "workdir": args.workdir,
                    "fasta_path": args.fasta,
                    "workflow_path": args.workflow,
                    "thread_num": args.threads,
                }
            )
            common = {
                "run_id": args.run_id,
                "fdr_threshold": args.global_fdr,
            }

            if args.data_type == 'wiff2mzml':
                if not args.wiff_mzml:
                    parser.error("fp-msdt --instrument wiff requires --wiff-mzml")
                gen_wiff_fragpipe_msdt(
                    rawspectrum_parquet,
                    args.file,
                    pin_file,
                    target_tsv,
                    decoy_tsv,
                    final_output_file,
                    True,
                    **common,
                )
            else:
                gen_mzml_fragpipe_msdt(
                    rawspectrum_parquet,
                    pin_file,
                    final_output_file,
                    True,
                    target_tsv,
                    decoy_tsv,
                    **common,
                )

            return 0
        if args.command == "fp-msdt":
            common = {
                "run_id": args.run_id,
                "fdr_threshold": args.global_fdr,
            }
            if args.instrument == "wiff":
                if not args.wiff_mzml:
                    parser.error("fp-msdt --instrument wiff requires --wiff-mzml")
                state = gen_wiff_fragpipe_msdt(
                    args.raw_spectrum,
                    args.wiff_mzml,
                    args.pin,
                    args.target_tsv,
                    args.decoy_tsv,
                    args.output,
                    not args.no_unify_residue,
                    **common,
                )
            else:
                state = gen_mzml_fragpipe_msdt(
                    args.raw_spectrum,
                    args.pin,
                    args.output,
                    not args.no_unify_residue,
                    args.target_tsv,
                    args.decoy_tsv,
                    **common,
                )
            return 0 if state in (0, 1) else 1
        if args.command == "global-percolator":
            run_global_percolator(
                _pin_mapping(args.pin),
                args.percolator_exe,
                args.output_dir,
                threads=args.threads,
            )
            return 0
        if args.legacy_config:
            return run_config(args.legacy_config)
        parser.print_help()
        return 0
    except Exception as error:
        logger.error("MSDT-Converter failed: %s", error)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
