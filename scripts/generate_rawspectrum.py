import os
import logging
import subprocess
import pandas as pd

# Configure logger (accessible by caller)
logger = logging.getLogger(__name__)

deal_mzml_rawspectrum = "./linux_mzml_rawspectrum"
deal_tims_rawspectrum = "./linux_d_rawspectrum"
deal_wiff_rawspectrum = "./wiff_mzml_rawspecturm"


def generate_rawspectrum_fn(param):
    """
    Generate rawspectrum file
    Return values:
        0: Successfully generated
        1: Already exists, no need to generate
        2: Input file does not exist
        -1: Generation failed
    """
    data_type = param.get('data_type')
    raw_data_path = param.get('input')
    output_path = param.get('output')

    if not os.path.exists(raw_data_path):
        logger.error(f"Input file does not exist: {raw_data_path}")
        return 2

    if os.path.exists(output_path):
        logger.info(f"Output already exists, skipping: {output_path}")
        return 1
    
    # Create the output directory if it does not exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    if data_type == 'mzml':
        if not raw_data_path.endswith('.mzML'):
            logger.error(f"Input file extension is not .mzML: {raw_data_path}")
            return 2
        temp_output_path = output_path.replace('.parquet', '_temp.tsv')
        cmd = [deal_mzml_rawspectrum, raw_data_path, temp_output_path]

    elif data_type == 'tims':
        if not raw_data_path.endswith('.d'):
            logger.error(f"Input file extension is not .tims: {raw_data_path}")
            return 2
        temp_output_path = output_path.replace('.parquet', '_temp.tsv')
        cmd = [deal_tims_rawspectrum, raw_data_path, temp_output_path]

    elif data_type == 'wiff2mzml':
        if not raw_data_path.endswith('.mzML'):
            logger.error(f"Input file extension is not .mzML: {raw_data_path}")
            return 2
        temp_output_path = output_path.replace('.parquet', '_temp.tsv')
        cmd = [deal_wiff_rawspectrum, raw_data_path, temp_output_path]
        
    else:
        logger.error(f"Unknown data type: {data_type}")
        return -1

    try:
        logger.info(f"Generating: {output_path}")
        logger.info(f"Executing command: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        if data_type == 'tims':
            temp_df = pd.read_csv(temp_output_path, sep='\t', usecols=['scan','precursor_mz','rt','ion_mobility','mz_array','intensity_array'])
        else:
            temp_df = pd.read_csv(temp_output_path, sep='\t', usecols=['scan','precursor_mz','rt','mz_array','intensity_array'])
        temp_df.to_parquet(output_path)
        os.remove(temp_output_path)
        
        logger.info(f"Successfully generated: {output_path}")
        return 0
    except subprocess.CalledProcessError as e:
        logger.error(f"Generation failed: {e.stderr}")
        return -1
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        return -1