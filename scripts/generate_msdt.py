import os
import logging
import subprocess
import pandas as pd
import numpy as np
import re
from pathlib import Path
from typing import Callable

from scripts.percolator import (
    enrich_dataframe_with_percolator,
    normalize_modified_sequence,
    parse_psm_id,
)

# Configure logger
logger = logging.getLogger(__name__)

deal_mzml_rawspectrum = "./linux_mzml_rawspectrum"
deal_tims_rawspectrum = "./linux_d_rawspectrum"
deal_wiff_rawspectrum = "./wiff_mzml_rawspecturm"

residues_sage = {
    'C[+57.0216]': 'C[57.02]',
    'M[+15.9949]': 'M[15.99]',
    '[+42]-': 'n[42]',
    'N[+0.98]':'N[.98]',
    'Q[+0.98]':'Q[.98]'
}

residues_frag = {
    'M[147]': 'M[15.99]',
    'M[15.9949]': 'M[15.99]',
    'C[57.0215]': 'C[57.02]',
    'n[42.0106]': 'n[42]',
    'n[420106]': 'n[42]',
    'M[159949]': 'M[15.99]',
    'C[57.02][570215]': 'C[57.02]',
    'C[570215]': 'C[57.02]',
    'N[0.9800]': 'N[.98]',
    'Q[0.9800]': 'Q[.98]'
}

FRAGPIPE_PIN_REQUIRED_COLUMNS = [
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
    "Peptide",
    "Proteins",
]
FRAGPIPE_PIN_OPTIONAL_COLUMNS = [
    "unweighted_spectral_entropy",
    "delta_RT_loess",
]


def _read_fragpipe_pin(path):
    """Read a FragPipe PIN while tolerating optional MSBooster features."""
    frame = pd.read_csv(path, sep="\t")
    missing_required = sorted(
        set(FRAGPIPE_PIN_REQUIRED_COLUMNS).difference(frame.columns)
    )
    if missing_required:
        raise ValueError(
            f"FragPipe PIN {path} is missing required columns: "
            + ", ".join(missing_required)
        )
    missing_optional = [
        column
        for column in FRAGPIPE_PIN_OPTIONAL_COLUMNS
        if column not in frame.columns
    ]
    for column in missing_optional:
        frame[column] = np.nan
    if missing_optional:
        logger.warning(
            "FragPipe PIN %s is missing optional MSBooster columns: %s",
            path,
            ", ".join(missing_optional),
        )
    return frame[FRAGPIPE_PIN_REQUIRED_COLUMNS + FRAGPIPE_PIN_OPTIONAL_COLUMNS]

def keep_uppercase(s: str) -> str:
    """Remove all characters in the string that are not capital letters."""
    return re.sub(r'[^A-Z]', '', s)

def clean_psm_func(peptide, residues_dict):
    for key, value in residues_dict.items():
        if value not in peptide:
            peptide = peptide.replace(key, value)
    return peptide
    
def change_wiff_scan(right_scan_path, wrong_scan_path):
    right_df = pd.read_parquet(right_scan_path)
    wrong_df = pd.read_csv(wrong_scan_path, sep='\t', usecols=['scan'])
    right_df['scan_sr'] = list(wrong_df['scan'])
    return right_df
 
def gen_mzml_tims_sage_msdt(raw_data_path, search_result_path, output_path, unify_residue):
    try:
        raw_df = pd.read_parquet(raw_data_path)
        raw_df = raw_df[['scan','precursor_mz','rt','mz_array','intensity_array']]
        raw_df = raw_df.dropna(subset=['scan', 'mz_array','intensity_array'])
        raw_df['scan'] = raw_df['scan'].astype(int)
        raw_df['mz_array'] = raw_df['mz_array'].str.split(',').map(lambda x: np.array(x, dtype='float32'))
        raw_df['intensity_array'] = raw_df['intensity_array'].str.split(',').map(lambda x: np.array(x, dtype='float32'))
        
        sage_df = pd.read_csv(search_result_path,sep='\t',usecols=['peptide','scannr','label','matched_peaks','peptide_q','protein_q','charge','predicted_rt','ion_mobility','delta_rt_model','sage_discriminant_score','spectrum_q','proteins'])
        sage_df = sage_df.dropna(subset=['scannr'])
        try:
            sage_df['scannr'] = sage_df['scannr'].astype(int)
        except:
            sage_df['scannr'] = sage_df['scannr'].apply(lambda x: x.split('=')[-1]).astype(int)
        decoy_df = sage_df[sage_df['label']==-1]
        target_df = sage_df[sage_df['label']==1]
        identify_target_df = target_df[target_df['spectrum_q']<=0.01]
        sage_df_need = pd.concat([identify_target_df, decoy_df], axis=0, ignore_index=True)
        sage_df_need = sage_df_need.rename(columns={'scannr':'scan','delta_rt_model':'delta_rt'})
        if unify_residue:
            sage_df_need['precursor_sequence'] = sage_df_need['peptide'].apply(lambda x: clean_psm_func(x,residues_sage))
        else:
            sage_df_need['precursor_sequence'] = sage_df_need['peptide']
        sage_df_need['cleaned_sequence'] = sage_df_need['precursor_sequence'].apply(lambda x: keep_uppercase(x))
        # sage_df_need['cleaned_sequence'] = sage_df_need['precursor_sequence'].str.replace('n[42]', '').str.replace('N[.98]', 'N').str.replace('Q[.98]', 'Q').str.replace('M[15.99]', 'M').str.replace('C[57.02]', 'C')
        sage_df_need['sequence_len'] = sage_df_need['cleaned_sequence'].apply(len)
        sage_df_need = sage_df_need[(sage_df_need['sequence_len']<=50)&(sage_df_need['sequence_len']>=7)]
        sage_df_need = sage_df_need[(sage_df_need['charge']<=5)&(sage_df_need['charge']>=2)]
        sage_df_need['label'] = sage_df_need['label'].replace(-1, 0)
        sage_df_need['ion_mobility'] = sage_df_need['ion_mobility'].fillna(0)
        sage_df_need = sage_df_need[['scan','precursor_sequence','proteins','label','matched_peaks','peptide_q','protein_q','charge','predicted_rt','ion_mobility','delta_rt','spectrum_q','sage_discriminant_score']]
        resultdf_grouped = sage_df_need.groupby('scan').agg({
                'precursor_sequence': lambda x: list(x),
                'proteins': lambda x: list(x),
                'label': lambda x: list(x),
                'charge': lambda x: list(x),
                'matched_peaks': lambda x: list(x),
                'peptide_q': lambda x: list(x),
                'protein_q': lambda x: list(x),
                'predicted_rt': lambda x: list(x),
                'ion_mobility': lambda x: list(x),
                'delta_rt': lambda x: list(x),
                'spectrum_q': lambda x: list(x),
                'sage_discriminant_score': lambda x: list(x)
            }).reset_index()
        resultdf_grouped['spectrum_q'] = resultdf_grouped['spectrum_q'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['sage_discriminant_score'] = resultdf_grouped['sage_discriminant_score'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['label'] = resultdf_grouped['label'].apply(lambda x: [np.int8(i) for i in x])
        resultdf_grouped['charge'] = resultdf_grouped['charge'].apply(lambda x: [np.int8(i) for i in x])
        resultdf_grouped['matched_peaks'] = resultdf_grouped['matched_peaks'].apply(lambda x: [np.int32(i) for i in x])
        resultdf_grouped['peptide_q'] = resultdf_grouped['peptide_q'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['protein_q'] = resultdf_grouped['protein_q'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['predicted_rt'] = resultdf_grouped['predicted_rt'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['ion_mobility'] = resultdf_grouped['ion_mobility'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['delta_rt'] = resultdf_grouped['delta_rt'].apply(lambda x: [np.float32(i) for i in x])

        parquet_df = resultdf_grouped.merge(raw_df, on='scan', how='inner')
        assert len(parquet_df) == len(resultdf_grouped)

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        parquet_df.to_parquet(output_path)
        return 0
    except Exception as e:
        logger.error(f"Error occurs when generate {output_path}: {e}")
        return -1
    
def gen_mzml_fragpipe_msdt(
    raw_data_path,
    fp_pin_path,
    output_path,
    unify_residue,
    percolator_target_path=None,
    percolator_decoy_path=None,
    run_id=None,
    fdr_threshold=None,
):
    try:
        raw_df = pd.read_parquet(raw_data_path)
        if 'ion_mobility' in raw_df.columns:
            raw_df = raw_df[['scan','precursor_mz','rt','ion_mobility','mz_array','intensity_array']]
        else:
            raw_df = raw_df[['scan','precursor_mz','rt','mz_array','intensity_array']]
        raw_df = raw_df.dropna(subset=['scan', 'mz_array','intensity_array'])
        raw_df['scan'] = raw_df['scan'].astype(int)
        raw_df['mz_array'] = raw_df['mz_array'].map(lambda x: np.asarray(x.split(',') if isinstance(x, str) else x, dtype='float32'))
        raw_df['intensity_array'] = raw_df['intensity_array'].map(lambda x: np.asarray(x.split(',') if isinstance(x, str) else x, dtype='float32'))

        # read fp_sr decoy
        fp_sr_df = _read_fragpipe_pin(fp_pin_path)
        fp_sr_df = fp_sr_df.rename(columns={'ScanNr': 'scan', 'Label':'label', 'Proteins': 'proteins'})
        fp_sr_df['scan'] = fp_sr_df['scan'].astype(int)
        fp_sr_df['label'] = fp_sr_df['label'].replace(-1, 0).astype('int8')
        fp_sr_df['charge'] = fp_sr_df['SpecId'].map(lambda value: parse_psm_id(value)[1])
        
        if unify_residue:
            fp_sr_df['precursor_sequence'] = fp_sr_df['Peptide'].map(normalize_modified_sequence)
        else:
            fp_sr_df['precursor_sequence'] = fp_sr_df['Peptide'].map(lambda value: re.sub(r'^[A-Za-z-]\.(.+)\.[A-Za-z-]$', r'\1', str(value)))
        fp_sr_df = fp_sr_df[['scan', 'label', 'charge', 'ExpMass', 'retentiontime', 'rank', 'isotope_errors', 'hyperscore', 'delta_hyperscore', 'matched_ion_num', 'ion_series', 'unweighted_spectral_entropy', 'delta_RT_loess', 'precursor_sequence', 'proteins']]

        fp_parquet_df = fp_sr_df.merge(raw_df, on='scan', how='left', validate='many_to_one', indicator=True)
        unmatched = fp_parquet_df['_merge'] != 'both'
        if unmatched.any():
            examples = fp_parquet_df.loc[unmatched, 'scan'].head(5).tolist()
            raise ValueError(f"FragPipe scans missing from raw spectra: {examples}")
        fp_parquet_df = fp_parquet_df.drop(columns='_merge')
        if percolator_target_path and percolator_decoy_path:
            fp_parquet_df, _ = enrich_dataframe_with_percolator(
                fp_parquet_df,
                percolator_target_path,
                percolator_decoy_path,
                run_id=run_id,
                fdr_threshold=fdr_threshold,
            )
        elif percolator_target_path or percolator_decoy_path:
            raise ValueError("Both Percolator target and decoy TSV paths are required")
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fp_parquet_df.to_parquet(output_path, index=False)
        return 0
    except Exception as e:
        logger.error(f"Error occurs when generate {output_path}: {e}")
        return -1


def gen_wiff_fragpipe_msdt(
    raw_data_path,
    wiff_mzml_path,
    fp_pin_path,
    percolator_target_path,
    percolator_decoy_path,
    output_path,
    unify_residue=True,
    *,
    run_id=None,
    fdr_threshold=None,
    mzml_extractor=deal_mzml_rawspectrum,
    runner: Callable = subprocess.run,
):
    """Generate a Percolator-enriched FP MSDT while preserving WIFF scans."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    search_scan_tsv = output.parent / f"{Path(wiff_mzml_path).stem}_search_scans.tsv"
    try:
        runner(
            [str(mzml_extractor), str(wiff_mzml_path), str(search_scan_tsv)],
            check=True,
            capture_output=True,
            text=True,
        )
        if not search_scan_tsv.is_file():
            raise RuntimeError(
                f"mzML extractor did not create search-scan TSV: {search_scan_tsv}"
            )

        raw_df = pd.read_parquet(raw_data_path)
        required_raw = {
            "scan",
            "precursor_mz",
            "rt",
            "mz_array",
            "intensity_array",
        }
        missing_raw = sorted(required_raw.difference(raw_df.columns))
        if missing_raw:
            raise ValueError(
                "WIFF raw-spectrum Parquet is missing columns: "
                + ", ".join(missing_raw)
            )
        search_scans = pd.read_csv(search_scan_tsv, sep="\t")
        required_search = {"scan", "precursor_mz", "rt"}
        missing_search = sorted(required_search.difference(search_scans.columns))
        if missing_search:
            raise ValueError(
                "mzML search-scan TSV is missing validation columns: "
                + ", ".join(missing_search)
            )
        search_scans = search_scans[["scan", "precursor_mz", "rt"]].copy()
        if len(raw_df) != len(search_scans):
            raise ValueError(
                "WIFF native/search scan row counts differ: "
                f"{len(raw_df)} != {len(search_scans)}"
            )
        if raw_df["scan"].duplicated().any():
            raise ValueError("WIFF native scan values are not unique")
        if search_scans["scan"].duplicated().any():
            raise ValueError("WIFF search scan values are not unique")

        native_mz = pd.to_numeric(raw_df["precursor_mz"], errors="raise").to_numpy()
        search_mz = pd.to_numeric(
            search_scans["precursor_mz"], errors="raise"
        ).to_numpy()
        native_rt = pd.to_numeric(raw_df["rt"], errors="raise").to_numpy()
        search_rt = pd.to_numeric(search_scans["rt"], errors="raise").to_numpy()
        invalid_spectrum_values = ~(
            np.isfinite(native_mz)
            & np.isfinite(search_mz)
            & np.isfinite(native_rt)
            & np.isfinite(search_rt)
        )
        spectrum_mismatch = invalid_spectrum_values | ~(
            np.isclose(native_mz, search_mz, rtol=0, atol=0.01)
            & np.isclose(native_rt, search_rt, rtol=0, atol=0.05)
        )
        if spectrum_mismatch.any():
            mismatch_positions = np.flatnonzero(spectrum_mismatch)[:5]
            examples = [
                f"row {position}: native scan {raw_df.iloc[position]['scan']} "
                f"vs search scan {search_scans.iloc[position]['scan']}"
                for position in mismatch_positions
            ]
            raise ValueError(
                "WIFF native/search spectrum order does not match within "
                "precursor_mz=0.01 and rt=0.05 tolerances: "
                + "; ".join(examples)
            )

        raw_df = raw_df[
            ["scan", "precursor_mz", "rt", "mz_array", "intensity_array"]
        ].copy()
        native_spectrum_index = pd.to_numeric(raw_df["scan"], errors="raise")
        if not np.equal(native_spectrum_index % 1, 0).all():
            raise ValueError("WIFF raw-spectrum scan values must be integer mzML indices")
        # MSFragger reports the mzML spectrum index as a 1-based ScanNr. SCIEX
        # native IDs (sample/period/cycle/experiment) are deliberately used only
        # for order validation because cycle alone is not unique.
        raw_df["search_scan"] = native_spectrum_index.astype(int) + 1
        raw_df = raw_df.dropna(subset=["scan", "mz_array", "intensity_array"])
        raw_df["scan"] = raw_df["scan"].astype(int)

        def to_float_array(value):
            if isinstance(value, str):
                value = value.split(",")
            return np.asarray(value, dtype="float32")

        raw_df["mz_array"] = raw_df["mz_array"].map(to_float_array)
        raw_df["intensity_array"] = raw_df["intensity_array"].map(
            to_float_array
        )

        fp_df = _read_fragpipe_pin(fp_pin_path)
        fp_df = fp_df.rename(
            columns={
                "ScanNr": "search_scan",
                "Label": "label",
                "Proteins": "proteins",
            }
        )
        fp_df["search_scan"] = fp_df["search_scan"].astype(int)
        fp_df["label"] = fp_df["label"].replace(-1, 0).astype("int8")
        fp_df["charge"] = fp_df["SpecId"].map(
            lambda value: parse_psm_id(value)[1]
        )
        if unify_residue:
            fp_df["precursor_sequence"] = fp_df["Peptide"].map(
                normalize_modified_sequence
            )
        else:
            fp_df["precursor_sequence"] = fp_df["Peptide"].map(
                lambda value: re.sub(
                    r"^[A-Za-z-]\.(.+)\.[A-Za-z-]$", r"\1", str(value)
                )
            )
        fp_df = fp_df[
            [
                "search_scan",
                "label",
                "charge",
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
                "precursor_sequence",
                "proteins",
            ]
        ]
        base_msdt = fp_df.merge(
            raw_df,
            on="search_scan",
            how="left",
            validate="many_to_one",
            indicator=True,
        )
        unmatched_scan = base_msdt["_merge"] != "both"
        if unmatched_scan.any():
            examples = (
                base_msdt.loc[unmatched_scan, "search_scan"].head(5).tolist()
            )
            raise ValueError(
                "FragPipe search scans missing from WIFF scan map: "
                + ", ".join(str(value) for value in examples)
            )
        base_msdt = base_msdt.drop(columns="_merge")
        enriched, _ = enrich_dataframe_with_percolator(
            base_msdt,
            percolator_target_path,
            percolator_decoy_path,
            run_id=run_id,
            scan_column="search_scan",
            fdr_threshold=fdr_threshold,
        )
        enriched = enriched.drop(columns="search_scan")
        enriched.to_parquet(output, index=False)
        return 0
    except Exception as error:
        logger.exception(f"Error occurs when generating WIFF FP MSDT {output}: {error}")
        return -1
    finally:
        search_scan_tsv.unlink(missing_ok=True)
    
def gen_wiff_sage_msdt(raw_data_path, wiff_mzml_path, search_result_path, output_path, unify_residue):
    try:
        wiff_fn = wiff_mzml_path.split('/')[-1][:-5]
        wrong_temp_raw = f"{'/'.join(output_path.split('/')[:-1])}/{wiff_fn}_temp_raw.tsv"
        cmd = [deal_mzml_rawspectrum, wiff_mzml_path, wrong_temp_raw]
        subprocess.run(cmd, capture_output=True, text=True)
        
        if os.path.exists(wrong_temp_raw):
            raw_df = change_wiff_scan(raw_data_path, wrong_temp_raw)
            raw_df = raw_df.dropna(subset=['scan', 'mz_array','intensity_array'])
            raw_df['scan'] = raw_df['scan'].astype(int)
            raw_df['mz_array'] = raw_df['mz_array'].str.split(',').map(lambda x: np.array(x, dtype='float32'))
            raw_df['intensity_array'] = raw_df['intensity_array'].str.split(',').map(lambda x: np.array(x, dtype='float32'))
        else:
            logger.error(f"fail to generate wiff temp rawspectrum: {wrong_temp_raw}")
            return -1

        sage_df = pd.read_csv(search_result_path,sep='\t',usecols=['peptide','scannr','label','matched_peaks','peptide_q','protein_q','charge','predicted_rt','ion_mobility','delta_rt_model','sage_discriminant_score','spectrum_q','proteins'])
        sage_df = sage_df.dropna(subset=['scannr'])
        decoy_df = sage_df[sage_df['label']==-1]
        target_df = sage_df[sage_df['label']==1]
        identify_target_df = target_df[target_df['spectrum_q']<=0.01]
        sage_df_need = pd.concat([identify_target_df, decoy_df], axis=0, ignore_index=True)

        sage_df_need = sage_df_need.rename(columns={'scannr':'scan_sr','delta_rt_model':'delta_rt'})
        if unify_residue:
            sage_df_need['precursor_sequence'] = sage_df_need['peptide'].apply(lambda x: clean_psm_func(x,residues_sage))
        else:
            sage_df_need['precursor_sequence'] = sage_df_need['peptide']
        sage_df_need['cleaned_sequence'] = sage_df_need['precursor_sequence'].apply(lambda x: keep_uppercase(x))
        sage_df_need['sequence_len'] = sage_df_need['cleaned_sequence'].apply(len)
        sage_df_need = sage_df_need[(sage_df_need['sequence_len']<=50)&(sage_df_need['sequence_len']>=7)]
        sage_df_need = sage_df_need[(sage_df_need['charge']<=5)&(sage_df_need['charge']>=2)]
        sage_df_need['label'] = sage_df_need['label'].replace(-1, 0)
        sage_df_need['ion_mobility'] = sage_df_need['ion_mobility'].fillna(0)
        sage_df_need = sage_df_need[['scan_sr','precursor_sequence','proteins','label','matched_peaks','peptide_q','protein_q','charge','predicted_rt','ion_mobility','delta_rt','spectrum_q','sage_discriminant_score']]
        resultdf_grouped = sage_df_need.groupby('scan_sr').agg({
                'precursor_sequence': lambda x: list(x),
                'proteins': lambda x: list(x),
                'label': lambda x: list(x),
                'charge': lambda x: list(x),
                'matched_peaks': lambda x: list(x),
                'peptide_q': lambda x: list(x),
                'protein_q': lambda x: list(x),
                'predicted_rt': lambda x: list(x),
                'ion_mobility': lambda x: list(x),
                'delta_rt': lambda x: list(x),
                'spectrum_q': lambda x: list(x),
                'sage_discriminant_score': lambda x: list(x)
            }).reset_index()
        resultdf_grouped['spectrum_q'] = resultdf_grouped['spectrum_q'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['sage_discriminant_score'] = resultdf_grouped['sage_discriminant_score'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['label'] = resultdf_grouped['label'].apply(lambda x: [np.int8(i) for i in x])
        resultdf_grouped['charge'] = resultdf_grouped['charge'].apply(lambda x: [np.int8(i) for i in x])
        resultdf_grouped['matched_peaks'] = resultdf_grouped['matched_peaks'].apply(lambda x: [np.int32(i) for i in x])
        resultdf_grouped['peptide_q'] = resultdf_grouped['peptide_q'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['protein_q'] = resultdf_grouped['protein_q'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['predicted_rt'] = resultdf_grouped['predicted_rt'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['ion_mobility'] = resultdf_grouped['ion_mobility'].apply(lambda x: [np.float32(i) for i in x])
        resultdf_grouped['delta_rt'] = resultdf_grouped['delta_rt'].apply(lambda x: [np.float32(i) for i in x])

        parquet_df = resultdf_grouped.merge(raw_df, on='scan_sr', how='inner')
        parquet_df = parquet_df.drop('scan_sr', axis=1)
        assert len(parquet_df) == len(resultdf_grouped)

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        parquet_df.to_parquet(output_path)
        return 0
    except Exception as e:
        logger.error(f"Error occurs when generate {output_path}: {e}")
        return -1
    
def generate_msdt_fn(param):
    """Generate all enabled MSDT variants, including WIFF FragPipe output."""
    states = []
    outputs = []

    def missing(paths):
        absent = [str(path) for path in paths if not path or not Path(path).exists()]
        if absent:
            logger.error("Missing MSDT inputs: %s", ", ".join(absent))
            states.append(2)
            return True
        return False

    tims = param.get("tims", {})
    if tims.get("need_tims", False):
        output = tims.get("output")
        outputs.append(output)
        if output and Path(output).exists():
            states.append(1)
        elif not missing(
            [tims.get("rawspectrum_path"), tims.get("sage_search_result_path")]
        ):
            states.append(
                gen_mzml_tims_sage_msdt(
                    tims["rawspectrum_path"],
                    tims["sage_search_result_path"],
                    output,
                    tims.get("unify_residue", True),
                )
            )

    mzml = param.get("mzml", {})
    if mzml.get("need_mzml", False):
        raw_path = mzml.get("rawspectrum_path")
        if mzml.get("need_sage", False):
            output = mzml.get("sage_output")
            outputs.append(output)
            if output and Path(output).exists():
                states.append(1)
            elif not missing([raw_path, mzml.get("sage_search_result_path")]):
                states.append(
                    gen_mzml_tims_sage_msdt(
                        raw_path,
                        mzml["sage_search_result_path"],
                        output,
                        mzml.get("sage_unify_residue", True),
                    )
                )
        if mzml.get("need_fragpipe", False):
            output = mzml.get("fp_output")
            outputs.append(output)
            percolator_target = mzml.get("percolator_target_path")
            percolator_decoy = mzml.get("percolator_decoy_path")
            required = [raw_path, mzml.get("fp_pin_path")]
            if percolator_target or percolator_decoy:
                required.extend([percolator_target, percolator_decoy])
            if output and Path(output).exists():
                states.append(1)
            elif not missing(required):
                states.append(
                    gen_mzml_fragpipe_msdt(
                        raw_path,
                        mzml["fp_pin_path"],
                        output,
                        mzml.get("fp_unify_residue", True),
                        percolator_target,
                        percolator_decoy,
                        mzml.get("run_id"),
                        mzml.get("fdr_threshold"),
                    )
                )

    wiff = param.get("wiff", {})
    if wiff.get("need_wiff", False):
        raw_path = wiff.get("rawspectrum_path")
        mzml_path = wiff.get("wiff_mzml_path")
        need_sage = wiff.get("need_sage", "need_fragpipe" not in wiff)
        if need_sage:
            output = wiff.get("sage_output", wiff.get("output"))
            outputs.append(output)
            if output and Path(output).exists():
                states.append(1)
            elif not missing(
                [raw_path, mzml_path, wiff.get("sage_search_result_path")]
            ):
                states.append(
                    gen_wiff_sage_msdt(
                        raw_path,
                        mzml_path,
                        wiff["sage_search_result_path"],
                        output,
                        wiff.get(
                            "sage_unify_residue", wiff.get("unify_residue", True)
                        ),
                    )
                )
        if wiff.get("need_fragpipe", False):
            output = wiff.get("fp_output")
            outputs.append(output)
            required = [
                raw_path,
                mzml_path,
                wiff.get("fp_pin_path"),
                wiff.get("percolator_target_path"),
                wiff.get("percolator_decoy_path"),
            ]
            if output and Path(output).exists():
                states.append(1)
            elif not missing(required):
                states.append(
                    gen_wiff_fragpipe_msdt(
                        raw_path,
                        mzml_path,
                        wiff["fp_pin_path"],
                        wiff["percolator_target_path"],
                        wiff["percolator_decoy_path"],
                        output,
                        wiff.get("fp_unify_residue", True),
                        run_id=wiff.get("run_id"),
                        fdr_threshold=wiff.get("fdr_threshold"),
                    )
                )

    if not states:
        logger.info("No MSDT variants enabled")
        return 0
    if all(state in (0, 1) for state in states):
        for output in outputs:
            logger.info("generate MSDT success: %s", output)
        return 0
    logger.error("generate MSDT failed; states=%s", states)
    return -1
