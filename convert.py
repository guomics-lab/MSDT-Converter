import pandas as pd
import os
import numpy as np
import shutil
import argparse
import numpy as np
import re
import json
import logging

from scripts.generate_rawspectrum import generate_rawspectrum_fn
from scripts.generate_msdt import generate_msdt_fn
from scripts.mgf2parquet import mgf_to_parquet
from scripts.search_engine import generate_sage_search_result_fn, generate_fp_search_result_fn
from scripts.msdt2mgf import msdt2mgf


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# TODO:
# generate_rawspectrum:
#     1. mzml2rawspectrum
#     2. tims2rawspectrum
# generate_msdt:
#     1. sage + rawspectrum -> msdt
#     2. fp + rawspectrum -> msdt
# convert2msdt:
#     1. mgf -> msdt
#     2. mzIdentML -> msdt


parser = argparse.ArgumentParser()
parser.add_argument('-config', type=str, default="", help="convert config json")
args = parser.parse_args()

def load_config(config_path: str):
    """
    Read the configuration file and remove all fields starting with _comment
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path, "r", encoding="utf-8") as f:
        config = json.load(f)

    def remove_comments(obj):
        """Recursively remove fields that start with _comment"""
        if isinstance(obj, dict):
            return {k: remove_comments(v) for k, v in obj.items() if not k.startswith("_comment")}
        elif isinstance(obj, list):
            return [remove_comments(item) for item in obj]
        else:
            return obj

    return remove_comments(config)

def parse_config(cfg: dict):
    """
    Extract required parameters from the configuration dictionary.
    Returns a dict containing the steps to execute and related file paths.
    """
    steps = {}

    # Step 1: generate_rawspectrum
    if cfg.get("generate_rawspectrum", {}).get("need", False) == True:
        steps["generate_rawspectrum"] = {
            "data_type": cfg["generate_rawspectrum"]["data_type"],
            "input": cfg["generate_rawspectrum"]["data_path"],
            "output": cfg["generate_rawspectrum"]["output"]
        }

    # Step 2_1: generate sage search result
    if cfg.get("generate_sage_search_result", {}).get("need", False) == True:
        steps['generate_sage_search_result'] = {
            'workdir': cfg["generate_sage_search_result"]['workdir'],
            'fasta': cfg["generate_sage_search_result"]['fasta'],
            'data_path': cfg["generate_sage_search_result"]['data_path'],
            'config_path': cfg["generate_sage_search_result"]['config_path']
        }

    # Step 2_2: generate sage search result
    if cfg.get("generate_fragpipe_search_result", {}).get("need", False) == True:
        steps['generate_fragpipe_search_result'] = cfg.get("generate_fragpipe_search_result")

    # Step 3: generate_msdt
    if cfg.get("generate_msdt", {}).get("need", False) == True:
        steps["generate_msdt"] = {
            "tims": cfg["generate_msdt"]["tims"],
            "mzml": cfg["generate_msdt"]["mzml"],
            "wiff": cfg["generate_msdt"]["wiff"]
        }

    # Step 4: convert_2_msdt
    if cfg.get("convert_2_msdt", {}).get("mgf", {}).get("need", False) == True:
        steps["convert_2_msdt"] = {
            "mgf": cfg["convert_2_msdt"]["mgf"]
        }

    # Step 5: msdt2mgf
    if cfg.get("msdt_2_mgf", {}).get("need", False) == True:
        steps["msdt2mgf"] = {
            "msdt_path": cfg["msdt_2_mgf"]["msdt_path"],
            "output_path": cfg["msdt_2_mgf"]["output_path"]
        }

    return steps


if __name__ == "__main__":
    cfg = load_config(args.config)
    steps = parse_config(cfg)
    if steps == {}:
        logger.info("No steps to execute.")
    for step, params in steps.items():
        logger.info(f"Step to execute: {step}")
        for k, v in params.items():
            logger.info(f"    {k}: {v}")

    # call functions according to configured steps
    if "generate_rawspectrum" in steps:
        # do_generate_rawspectrum(steps["generate_rawspectrum"])
        logger.info("Calling generate_rawspectrum function")
        rawspectrum_state = generate_rawspectrum_fn(steps['generate_rawspectrum'])

    if "generate_sage_search_result" in steps:
        # do_generate_sage_search_result(steps["generate_sage_search_result"])
        logger.info("Calling generate_sage_search_result function")
        sage_state = generate_sage_search_result_fn(steps['generate_sage_search_result'])

    if "generate_fragpipe_search_result" in steps:
        sage_state = generate_fp_search_result_fn(steps['generate_fragpipe_search_result'])

    if "generate_msdt" in steps:
        # do_generate_msdt(steps["generate_msdt"])
        logger.info("Calling generate_msdt function")
        msdt_state = generate_msdt_fn(steps['generate_msdt'])

    if "convert_2_msdt" in steps:
        # do_convert_2_msdt(steps["convert_2_msdt"])
        logger.info("Calling convert_2_msdt function")
        convert_msdt_state = mgf_to_parquet(steps['convert_2_msdt']['mgf'])
    
    if "msdt2mgf" in steps:
        # do_msdt2mgf(steps["msdt2mgf"])
        logger.info("Calling msdt2mgf function")
        msdt2mgf_state = msdt2mgf(steps["msdt2mgf"])