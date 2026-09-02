"""Percolator PSM parsing and MSDT enrichment."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import subprocess
from typing import Any, Callable, Mapping, Sequence

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "PSMId",
    "score",
    "q-value",
    "posterior_error_prob",
    "peptide",
}

RESIDUE_REPLACEMENTS = {
    "M[147]": "M[15.99]",
    "M[15.9949]": "M[15.99]",
    "C[57.0215]": "C[57.02]",
    "n[42.0106]": "n[42]",
    "n[420106]": "n[42]",
    "M[159949]": "M[15.99]",
    "C[57.02][570215]": "C[57.02]",
    "C[570215]": "C[57.02]",
    "N[0.9800]": "N[.98]",
    "Q[0.9800]": "Q[.98]",
}


class DuplicatePsmIdError(ValueError):
    """Raised when a PSM key is not unique within one mass-spectrometry run."""


@dataclass(frozen=True)
class EnrichmentReport:
    """Auditable row counts from one Parquet enrichment operation."""

    input_rows: int
    percolator_rows: int
    matched_rows: int
    unmatched_input_rows: int
    output_rows: int


def normalize_modified_sequence(peptide: str) -> str:
    """Normalize a FragPipe/Percolator peptide to the MSDT residue format."""
    sequence = str(peptide).strip()
    flanked = re.fullmatch(r"[A-Za-z-]\.(.+)\.[A-Za-z-]", sequence)
    if flanked:
        sequence = flanked.group(1)
    for source in sorted(RESIDUE_REPLACEMENTS, key=len, reverse=True):
        sequence = sequence.replace(source, RESIDUE_REPLACEMENTS[source])
    return sequence


def parse_psm_id(psm_id: str) -> tuple[int, int]:
    """Extract scan and charge from a FragPipe Percolator PSMId."""
    identifier = str(psm_id).strip()
    match = re.search(
        r"\.(?P<scan>\d+)\.\d+\.(?P<charge>\d+)(?:_\d+)?$", identifier
    )
    if match is None:
        match = re.search(
            r"_(?P<scan>\d+)_(?P<charge>\d+)(?:_\d+)?$", identifier
        )
    if match is None:
        raise ValueError(
            f"Cannot parse scan and charge from Percolator PSMId: {identifier!r}"
        )
    return int(match.group("scan")), int(match.group("charge"))


def make_psm_id(scan: Any, charge: Any, modified_sequence: Any) -> str:
    """Build the v2 PSM key from scan, charge and modified sequence."""
    return f"{int(scan)}_{int(charge)}_{modified_sequence}"


def _read_psms(
    path: str | Path, label: int, run_id: str | None = None
) -> pd.DataFrame:
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"Percolator TSV not found: {source}")
    frame = pd.read_csv(source, sep="\t")
    missing = sorted(REQUIRED_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError(
            f"Percolator TSV {source} is missing columns: {', '.join(missing)}"
        )
    if run_id is not None:
        prefix = f"{run_id}::"
        frame = frame[frame["PSMId"].astype(str).str.startswith(prefix)].copy()

    scan_charge = frame["PSMId"].map(parse_psm_id)
    frame["scan"] = scan_charge.map(lambda value: value[0])
    frame["charge"] = scan_charge.map(lambda value: value[1])
    frame["modified_sequence"] = frame["peptide"].map(
        normalize_modified_sequence
    )
    frame["psm_id"] = [
        make_psm_id(scan, charge, sequence)
        for scan, charge, sequence in zip(
            frame["scan"], frame["charge"], frame["modified_sequence"]
        )
    ]
    frame["label"] = np.int8(label)
    frame["score"] = pd.to_numeric(frame["score"], errors="raise")
    frame["q-value"] = pd.to_numeric(frame["q-value"], errors="raise")
    frame["PEP"] = pd.to_numeric(
        frame["posterior_error_prob"], errors="raise"
    )
    metric_columns = ["score", "q-value", "PEP"]
    finite_metrics = np.isfinite(frame[metric_columns].to_numpy(dtype=float))
    if frame[metric_columns].isna().any().any() or not finite_metrics.all():
        invalid_columns = [
            column
            for column in metric_columns
            if frame[column].isna().any()
            or not np.isfinite(frame[column].to_numpy(dtype=float)).all()
        ]
        raise ValueError(
            f"Percolator TSV {source} contains empty or non-finite values in: "
            + ", ".join(invalid_columns)
        )
    return frame[
        [
            "psm_id",
            "scan",
            "charge",
            "modified_sequence",
            "label",
            "score",
            "q-value",
            "PEP",
        ]
    ]


def load_percolator_psms(
    target_path: str | Path,
    decoy_path: str | Path,
    *,
    run_id: str | None = None,
) -> pd.DataFrame:
    """Load target and decoy Percolator PSM tables into one canonical table."""
    psms = pd.concat(
        [
            _read_psms(target_path, 1, run_id),
            _read_psms(decoy_path, 0, run_id),
        ],
        ignore_index=True,
    )
    duplicate_mask = psms["psm_id"].duplicated(keep=False)
    if duplicate_mask.any():
        examples = (
            psms.loc[duplicate_mask, "psm_id"].drop_duplicates().head(5).tolist()
        )
        raise DuplicatePsmIdError(
            "Duplicate psm_id values found after merging target and decoy TSVs: "
            + ", ".join(examples)
        )
    return psms


def enrich_dataframe_with_percolator(
    parquet: pd.DataFrame,
    target_path: str | Path,
    decoy_path: str | Path,
    *,
    run_id: str | None = None,
    scan_column: str = "scan",
    fdr_threshold: float | None = None,
) -> tuple[pd.DataFrame, EnrichmentReport]:
    """Enrich an in-memory MSDT table using a selectable search-scan column."""
    required = {scan_column, "charge", "precursor_sequence", "label"}
    missing = sorted(required.difference(parquet.columns))
    if missing:
        raise ValueError(
            "MSDT Parquet is missing columns: " + ", ".join(missing)
        )

    parquet = parquet.drop(
        columns=[
            column
            for column in ("psm_id", "score", "q-value", "PEP")
            if column in parquet.columns
        ]
    ).copy()
    parquet["psm_id"] = [
        make_psm_id(scan, charge, sequence)
        for scan, charge, sequence in zip(
            parquet[scan_column],
            parquet["charge"],
            parquet["precursor_sequence"],
        )
    ]
    duplicate_mask = parquet["psm_id"].duplicated(keep=False)
    if duplicate_mask.any():
        examples = (
            parquet.loc[duplicate_mask, "psm_id"].drop_duplicates().head(5).tolist()
        )
        raise DuplicatePsmIdError(
            "Duplicate psm_id values found in MSDT Parquet: "
            + ", ".join(examples)
        )

    psms = load_percolator_psms(target_path, decoy_path, run_id=run_id)
    if fdr_threshold is not None:
        if not 0 <= fdr_threshold <= 1:
            raise ValueError("fdr_threshold must be between 0 and 1")

    parquet_ids = set(parquet["psm_id"])
    percolator_ids = set(psms["psm_id"])
    missing_from_parquet = psms.loc[
        ~psms["psm_id"].isin(parquet_ids), "psm_id"
    ].head(5)
    if not missing_from_parquet.empty:
        raise ValueError(
            "Percolator PSMs missing from MSDT Parquet: "
            + ", ".join(missing_from_parquet.tolist())
        )
    missing_from_percolator = parquet.loc[
        ~parquet["psm_id"].isin(percolator_ids), "psm_id"
    ].head(5)
    if not missing_from_percolator.empty:
        raise ValueError(
            "MSDT Parquet PSMs missing from Percolator TSVs: "
            + ", ".join(missing_from_percolator.tolist())
        )

    matched = parquet.merge(
        psms,
        on="psm_id",
        how="left",
        validate="one_to_one",
        suffixes=("", "_percolator"),
    )

    label_mismatch = matched["label"] != matched["label_percolator"]
    if label_mismatch.any():
        examples = matched.loc[label_mismatch, "psm_id"].head(5).tolist()
        raise ValueError(
            "Target/decoy label mismatch for PSMs: " + ", ".join(examples)
        )

    matched_rows = len(matched)
    if fdr_threshold is not None:
        matched = matched[
            (matched["label_percolator"] == 0)
            | (matched["q-value"] <= fdr_threshold)
        ].copy()

    output = matched.drop(
        columns=[
            "psm_id",
            "scan_percolator",
            "charge_percolator",
            "modified_sequence",
            "label_percolator",
        ]
    )
    report = EnrichmentReport(
        input_rows=len(parquet),
        percolator_rows=len(psms),
        matched_rows=matched_rows,
        unmatched_input_rows=0,
        output_rows=len(output),
    )
    return output, report
