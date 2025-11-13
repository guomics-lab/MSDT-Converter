import os
import logging
import subprocess
import numpy as np
import pandas as pd
import sys

# Configure logger (accessible by caller)
logger = logging.getLogger(__name__)

# mgf template
mgf_template = """BEGIN IONS
TITLE={title}
PEPMASS={pepmass}
CHARGE={charge}
SCANS={scan}
RTINSECONDS={rt}
SEQ={seq}
LABEL={label}
PROTEIN={protein}
{spectrum}
END IONS

"""

def msdt2mgf(param):
    """
    Generate rawspectrum file
    Return values:
        0: Successfully generated
        1: Already exists, no need to generate
        2: Input file does not exist
        -1: Generation failed
    """
    msdt_path = param.get('msdt_path')
    output_path = param.get('output_path')

    if not os.path.exists(msdt_path):
        logger.error(f"Input file does not exist: {msdt_path}")
        return 2

    if os.path.exists(output_path):
        logger.info(f"Output already exists, skipping: {output_path}")
        return 1
    
    # Create the output directory if it does not exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    try:
        spec_df = pd.read_parquet(msdt_path)
        spec_df = spec_df[['scan','precursor_mz','charge','rt','label','mz_array','intensity_array','precursor_sequence','proteins']]
        if isinstance(spec_df['label'][0], np.ndarray):
            spec_df = spec_df.explode(['precursor_sequence','proteins','charge','label'])
        spec_df = spec_df.dropna(subset=['mz_array','intensity_array'])

        mgf_lines = []
        for _, row in spec_df.iterrows():
            try:
                mgf_lines.append(mgf_template.format(
                    title=msdt_path,
                    pepmass=row['precursor_mz'],
                    charge=str(int(row['charge'])) + '+',
                    scan=row['scan'],
                    rt=row['rt'] * 60,  # minutes to seconds
                    seq=row['precursor_sequence'],
                    label=row['label'],
                    protein=row['proteins'],
                    spectrum="\n".join(
                        f"{mz:.4f} {intensity:.2f}"  # control output precision
                        for mz, intensity in zip(
                            row['mz_array'], 
                            row['intensity_array']
                        )
                    )
                ))
            except Exception as e:
                print(f"Error processing row: {e}")
                sys.exit(1)
        
        with open(output_path, 'w') as f:
            f.writelines(mgf_lines)
        logger.info(f"Successfully generated: {output_path}")
        return 0

    except subprocess.CalledProcessError as e:
        logger.error(f"Generation failed: {e.stderr}")
        return -1
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        return -1
