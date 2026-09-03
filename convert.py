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
from scripts.search_engine import (
    generate_fp_search_result_fn,
    generate_sage_search_result_fn
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

    if "generate_msdt" in steps:
        msdt_params = copy.deepcopy(steps["generate_msdt"])
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
) -> int:
    """Run a JSON workflow with optional command-line overrides."""
    config = load_config(config_path)
    fragpipe = config.get("generate_fragpipe_search_result", {})
    if file_list is not None:
        fragpipe["file_list"] = file_list
        fragpipe.pop("data_path", None)
    if threads is not None:
        fragpipe["thread_num"] = threads
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


def _enrich_one(each_file_path, args):
    base_file_name = os.path.basename(each_file_path).removesuffix('.mzML')
    file_work_dir = os.path.join(args.workdir, base_file_name)
    os.makedirs(file_work_dir, exist_ok=True)
    rawspectrum_parquet = os.path.join(file_work_dir, base_file_name + "_rawspectrum.parquet")
    target_tsv = os.path.join(file_work_dir, 'exp', base_file_name + "_percolator_target_psms.tsv")
    decoy_tsv = os.path.join(file_work_dir, 'exp', base_file_name + "_percolator_decoy_psms.tsv")
    final_output_file = os.path.join(file_work_dir, base_file_name + "_fp_msdt_v2.parquet")
    pin_file = os.path.join(file_work_dir, 'exp', base_file_name + "_edited.pin")

    generate_rawspectrum_fn_param = {
        "data_type": args.data_type,
        "input": each_file_path,
        "output": rawspectrum_parquet,
    }
    generate_rawspectrum_fn(generate_rawspectrum_fn_param)
    generate_fp_search_result_fn(
        {
            "data_path": each_file_path,
            "workdir": file_work_dir,
            "fasta_path": args.fasta,
            "workflow_path": args.workflow,
            "thread_num": args.threads,
        }
    )
    common = {
    }
    if args.data_type == 'wiff2mzml':
        gen_wiff_fragpipe_msdt(
            rawspectrum_parquet,
            each_file_path,
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


def _read_file_list(txt_path):
    with open(txt_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    return [n.strip() for n in lines]

def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="MSDT-Converter v2")
    parser.add_argument('--version', action='version', version='v2.0')
    parser.add_argument(
        "-config", "--config", dest="legacy_config", help="legacy config JSON"
    )
    commands = parser.add_subparsers(dest="command")

    run = commands.add_parser("run", help="run a JSON pipeline")
    run.add_argument("--config", required=True)
    run.add_argument("--file-list")
    run.add_argument("--threads", type=int)

    enrich = commands.add_parser(
        "enrich", help="add Percolator fields to an FP MSDT Parquet"
    )
    enrich.add_argument("--file-list", required=True)
    enrich.add_argument("--data_type", required=True)
    enrich.add_argument("--workdir", required=True)
    enrich.add_argument("--fasta", required=True)
    enrich.add_argument("--workflow", required=True)
    enrich.add_argument("--threads", type=int, default=1)

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
            )

        if args.command == "enrich":
            file_list = _read_file_list(args.file_list)
            for abs_file_path in file_list:
                try:
                    _enrich_one(abs_file_path, args)
                except:
                    logger.exception(f"enrich file {abs_file_path} error")

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
