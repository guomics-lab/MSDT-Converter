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


@dataclass(frozen=True)
class GlobalPercolatorArtifacts:
    """Files produced by one pooled Percolator run."""

    combined_pin: Path
    target_tsv: Path
    decoy_tsv: Path


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


def enrich_parquet_with_percolator(
    parquet_path: str | Path,
    target_path: str | Path,
    decoy_path: str | Path,
    output_path: str | Path,
    *,
    run_id: str | None = None,
    fdr_threshold: float | None = None,
) -> EnrichmentReport:
    """Strictly match Percolator PSMs and add score, q-value and PEP."""
    parquet = pd.read_parquet(parquet_path)
    output, report = enrich_dataframe_with_percolator(
        parquet,
        target_path,
        decoy_path,
        run_id=run_id,
        fdr_threshold=fdr_threshold,
    )
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    output.to_parquet(destination, index=False)
    return report


def _combine_pin_files(pin_files: Mapping[str, str | Path], output: Path) -> None:
    if not pin_files:
        raise ValueError("At least one PIN file is required for global Percolator")

    expected_header = None
    default_direction = None
    combined_rows: list[str] = []
    next_global_scan = 1
    for run_id, pin_path in pin_files.items():
        if not run_id or "::" in run_id:
            raise ValueError(f"Invalid run_id for global Percolator: {run_id!r}")
        source = Path(pin_path)
        if not source.is_file():
            raise FileNotFoundError(f"PIN file not found: {source}")
        lines = source.read_text(encoding="utf-8").splitlines()
        if not lines:
            raise ValueError(f"PIN file is empty: {source}")
        if expected_header is None:
            expected_header = lines[0]
        elif lines[0] != expected_header:
            raise ValueError(
                f"PIN feature columns differ between files; incompatible file: {source}"
            )

        header_fields = lines[0].split("\t")
        try:
            spec_id_index = header_fields.index("SpecId")
            scan_index = header_fields.index("ScanNr")
        except ValueError as error:
            raise ValueError(
                f"PIN file must contain SpecId and ScanNr columns: {source}"
            ) from error
        run_scan_map: dict[int, int] = {}

        for line in lines[1:]:
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != len(header_fields):
                raise ValueError(f"Malformed PIN row in {source}: {line!r}")
            if fields[spec_id_index] == "DefaultDirection":
                if default_direction is None:
                    default_direction = line
                elif line != default_direction:
                    raise ValueError(
                        f"PIN DefaultDirection differs between files: {source}"
                    )
                continue
            try:
                source_scan = int(fields[scan_index])
            except ValueError as error:
                raise ValueError(
                    f"PIN ScanNr must be an integer in {source}: "
                    f"{fields[scan_index]!r}"
                ) from error
            if source_scan not in run_scan_map:
                run_scan_map[source_scan] = next_global_scan
                next_global_scan += 1
            fields[spec_id_index] = f"{run_id}::{fields[spec_id_index]}"
            fields[scan_index] = str(run_scan_map[source_scan])
            combined_rows.append("\t".join(fields))

    output_lines: list[str] = [expected_header or ""]
    if default_direction is not None:
        output_lines.append(default_direction)
    output_lines.extend(combined_rows)
    output.write_text("\n".join(output_lines) + "\n", encoding="utf-8")


def run_global_percolator(
    pin_files: Mapping[str, str | Path],
    percolator_executable: str | Path,
    output_dir: str | Path,
    *,
    threads: int = 1,
    runner: Callable[..., subprocess.CompletedProcess] = subprocess.run,
) -> GlobalPercolatorArtifacts:
    """Pool compatible PINs and run Percolator once for true global FDR."""
    executable = Path(percolator_executable)
    if not executable.is_file():
        raise FileNotFoundError(f"Percolator executable not found: {executable}")
    if threads < 1:
        raise ValueError("threads must be at least 1")

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    artifacts = GlobalPercolatorArtifacts(
        combined_pin=destination / "global_edited.pin",
        target_tsv=destination / "global_percolator_target_psms.tsv",
        decoy_tsv=destination / "global_percolator_decoy_psms.tsv",
    )
    _combine_pin_files(pin_files, artifacts.combined_pin)

    command: Sequence[str] = [
        str(executable),
        "--only-psms",
        "--no-terminate",
        "--post-processing-tdc",
        "--num-threads",
        str(threads),
        "--results-psms",
        str(artifacts.target_tsv),
        "--decoy-results-psms",
        str(artifacts.decoy_tsv),
        str(artifacts.combined_pin),
    ]
    runner(command, check=True, capture_output=True, text=True)
    missing = [
        path
        for path in (artifacts.target_tsv, artifacts.decoy_tsv)
        if not path.is_file()
    ]
    if missing:
        raise RuntimeError(
            "Percolator completed without expected TSV files: "
            + ", ".join(str(path) for path in missing)
        )
    return artifacts
