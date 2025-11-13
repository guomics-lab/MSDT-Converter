import os
import logging
import subprocess
import pandas as pd
import numpy as np
import re

# Configure logger
logger = logging.getLogger(__name__)

deal_mzml_rawspectrum = "../linux_mzml_rawspectrum"
deal_tims_rawspectrum = "../linux_d_rawspectrum"
deal_wiff_rawspectrum = "../wiff_mzml_rawspecturm"

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

def keep_uppercase(s: str) -> str:
    """Remove all characters in the string that are not capital letters."""
    return re.sub(r'[^A-Z]', '', s)

def clean_psm_func(peptide, residues_dict):
    for key, value in residues_dict.items():
        if value not in peptide:
            peptide = peptide.replace(key, value)
    return peptide
    
def change_wiff_scan(right_scan_path, wrong_scan_path):
    right_df = pd.read_csv(right_scan_path, sep='\t')
    wrong_df = pd.read_csv(wrong_scan_path, sep='\t', usecols=['scan','precursor_mz','rt'])
    wrong_df = wrong_df.rename(columns={'scan':'scan_sr'})
    merge_raw = right_df.merge(wrong_df, on=['precursor_mz','rt'], how='inner')
    assert len(merge_raw) == len(right_df), "ERROR: wiff rawspectrum match has different rows"
    return merge_raw
 
def gen_mzml_tims_sage_msdt(raw_data_path, search_result_path, output_path, unify_residue):
    try:
        raw_df = pd.read_csv(raw_data_path, sep='\t', usecols=['scan','precursor_mz','rt','mz_array','intensity_array'])
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
    
def gen_mzml_fragpipe_msdt(raw_data_path, fp_pin_path, output_path, unify_residue):
    try:
        raw_df = pd.read_csv(raw_data_path, sep='\t', usecols=['scan','precursor_mz','rt','mz_array','intensity_array'])
        raw_df = raw_df.dropna(subset=['scan', 'mz_array','intensity_array'])
        raw_df['scan'] = raw_df['scan'].astype(int)
        raw_df['mz_array'] = raw_df['mz_array'].str.split(',').map(lambda x: np.array(x, dtype='float32'))
        raw_df['intensity_array'] = raw_df['intensity_array'].str.split(',').map(lambda x: np.array(x, dtype='float32'))

        # read fp_sr decoy
        need_cols = ['SpecId', 'Label', 'ScanNr', 'ExpMass', 'retentiontime', 'rank', 'isotope_errors', 'hyperscore', 'delta_hyperscore', 'matched_ion_num', 'ion_series', 'unweighted_spectral_entropy', 'delta_RT_loess', 'Peptide', 'Proteins']
        fp_sr_df = pd.read_csv(fp_pin_path, sep='\t',usecols=need_cols)
        fp_sr_df = fp_sr_df.rename(columns={'ScanNr': 'scan', 'Label':'label', 'Proteins': 'proteins'})
        assert len(fp_sr_df) == len(set(fp_sr_df['scan'])), ""
        fp_sr_df['label'] = fp_sr_df['label'].replace(-1, 0)
        fp_sr_df['charge'] = fp_sr_df['SpecId'].apply(lambda x: int(x.split('.')[-1].split('_')[0]))
        
        def remove_trailing_numbers(s):
            return re.sub(r'\d+$', '', s)
        
        if unify_residue:
            fp_sr_df['precursor_sequence'] = fp_sr_df['Peptide'].apply(lambda x: clean_psm_func(x[2:-2], residues_frag))
        else:
            fp_sr_df['precursor_sequence'] = fp_sr_df['Peptide'].apply(lambda x: remove_trailing_numbers(x[2:-2]))
        fp_sr_df = fp_sr_df[['scan', 'label', 'charge', 'ExpMass', 'retentiontime', 'rank', 'isotope_errors', 'hyperscore', 'delta_hyperscore', 'matched_ion_num', 'ion_series', 'unweighted_spectral_entropy', 'delta_RT_loess', 'precursor_sequence', 'proteins']]

        fp_parquet_df = fp_sr_df.merge(raw_df, on='scan', how='inner')
        assert len(fp_parquet_df) == len(fp_sr_df)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        fp_parquet_df.to_parquet(output_path)
        return 0
    except Exception as e:
        logger.error(f"Error occurs when generate {output_path}: {e}")
        return -1
    
def gen_wiff_sage_msdt(raw_data_path, wiff_mzml_path, search_result_path, output_path, unify_residue):
    try:
        wiff_fn = wiff_mzml_path.split('/')[-1][:-5]
        wrong_temp_raw = f"{'/'.join(output_path.split('/')[:-1])}/{wiff_fn}_temp_raw.tsv"
        cmd = [deal_mzml_rawspectrum, wiff_mzml_path, wrong_temp_raw]
        subprocess.run(cmd, capture_output=True, text=True)
        
        if os.path.exists(wrong_temp_raw):
            merge_raw = change_wiff_scan(raw_data_path, wrong_temp_raw)
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

        parquet_df = resultdf_grouped.merge(merge_raw, on='scan_sr', how='inner')
        parquet_df = parquet_df.drop('scan_sr', axis=1)
        assert len(parquet_df) == len(resultdf_grouped)

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        parquet_df.to_parquet(output_path)
        return 0
    except Exception as e:
        logger.error(f"Error occurs when generate {output_path}: {e}")
        return -1
    
def generate_msdt_fn(param):
    """
    Generate a rawspectrum file.
    Return values:
        0: Successfully generated
        1: Already exists, no need to generate
        2: Input file does not exist
        -1: Generation failed
    """
    # tims
    tims_need = param.get('tims')['need_tims']
    tims_rawspectrum_path = param.get('tims')['rawspectrum_path']
    tims_sage_search_result_path = param.get('tims')['sage_search_result_path']
    tims_unify_residue = param.get('tims')['unify_residue']
    tims_output = param.get('tims')['output']
    tims_state = -1
    if tims_need:
        if os.path.exists(tims_output):
            logger.info(f"tims_output already done: {tims_output}, skip")
            tims_state = 1
        elif not os.path.exists(tims_sage_search_result_path):
            logger.error(f"miss tims_sage_search_result_path: {tims_sage_search_result_path}")
            tims_state = 2
        elif not os.path.exists(tims_rawspectrum_path):
            logger.error(f"miss tims_rawspectrum_path: {tims_rawspectrum_path}")
            tims_state = 2
        else:
            tims_state = gen_mzml_tims_sage_msdt(tims_rawspectrum_path, tims_sage_search_result_path, tims_output, tims_unify_residue)
        
    # mzml
    mzml_need = param.get('mzml')['need_mzml']
    mzml_need_sage = param.get('mzml')['need_sage']
    mzml_need_fragpipe = param.get('mzml')['need_fragpipe']
    mzml_rawspectrum_path = param.get('mzml')['rawspectrum_path']
    mzml_sage_search_result_path = param.get('mzml')['sage_search_result_path']
    mzml_fp_pin_path = param.get('mzml')['fp_pin_path']
    mzml_sage_unify_residue = param.get('mzml')['sage_unify_residue']
    mzml_fp_unify_residue = param.get('mzml')['fp_unify_residue']
    mzml_sage_output = param.get('mzml')['sage_output']
    mzml_fp_output = param.get('mzml')['fp_output']
    mzml_sage_state = -1
    mzml_fp_state = -1
    if mzml_need:
        if mzml_need_sage:
            if os.path.exists(mzml_sage_output):
                logger.info(f"mzml_sage_output already done: {mzml_sage_output}, skip")
                mzml_sage_state = 1
            elif not os.path.exists(mzml_rawspectrum_path):
                logger.error(f"miss mzml_rawspectrum_path: {mzml_rawspectrum_path}")
                mzml_sage_state = 2
            elif not os.path.exists(mzml_sage_search_result_path):
                logger.error(f"miss mzml_sage_search_result_path: {mzml_sage_search_result_path}")
                mzml_sage_state = 2
            else:
                mzml_sage_state = gen_mzml_tims_sage_msdt(mzml_rawspectrum_path, mzml_sage_search_result_path, mzml_sage_output, mzml_sage_unify_residue)
        else:
            mzml_sage_state = 0
                
        if mzml_need_fragpipe:
            
            if os.path.exists(mzml_fp_output):
                logger.info(f"mzml_fp_output already done: {mzml_fp_output}, skip")
                mzml_fp_state = 1
            elif not os.path.exists(mzml_rawspectrum_path):
                logger.error(f"miss mzml_rawspectrum_path: {mzml_rawspectrum_path}")
                mzml_fp_state = 2
            elif not os.path.exists(mzml_fp_pin_path):
                logger.error(f"miss mzml_fp_pin_path: {mzml_fp_pin_path}")
                mzml_fp_state = 2
            else:
                mzml_fp_state = gen_mzml_fragpipe_msdt(mzml_rawspectrum_path, mzml_fp_pin_path, mzml_fp_output, mzml_fp_unify_residue)
        else:
            mzml_fp_state = 0
    
    # wiff
    wiff_need = param.get('wiff')['need_wiff']
    wiff_mzml_path = param.get('wiff')['wiff_mzml_path']
    wiff_rawspectrum_path = param.get('wiff')['rawspectrum_path']
    wiff_sage_search_result_path = param.get('wiff')['sage_search_result_path']
    wiff_unify_residue = param.get('wiff')['unify_residue']
    wiff_output = param.get('wiff')['output']
    wiff_state = -1
    if wiff_need:
        if os.path.exists(wiff_output):
            logger.info(f"wiff_output already done: {wiff_output}, skip")
            wiff_state = 1
        elif not os.path.exists(wiff_rawspectrum_path):
            logger.error(f"miss wiff_rawspectrum_path: {wiff_rawspectrum_path}")
            wiff_state = 2
        elif not os.path.exists(wiff_mzml_path):
            logger.error(f"miss wiff_mzml_path: {wiff_mzml_path}")
            wiff_state = 2
        elif not os.path.exists(wiff_sage_search_result_path):
            logger.error(f"miss wiff_sage_search_result_path: {wiff_sage_search_result_path}")
            wiff_state = 2
        else:
            wiff_state = gen_wiff_sage_msdt(wiff_rawspectrum_path, wiff_mzml_path, wiff_sage_search_result_path, wiff_output, wiff_unify_residue)
    
    result_state = 0
    done_parquet_list = []
    if mzml_need:
        if mzml_sage_state != 0 or mzml_fp_state != 0:
            result_state = -1
        done_parquet_list.append(mzml_sage_output)
        done_parquet_list.append(mzml_fp_output)
    if tims_need:
        if tims_state != 0:
            result_state = -1
        done_parquet_list.append(tims_output)
    if wiff_need:
        if wiff_state != 0:
            result_state = -1
        done_parquet_list.append(wiff_output)
    
    if result_state == 0:
        logger.info(f"generate msdt success:")
        for i in done_parquet_list:
            logger.info(f"    generate msdt: {i} success")
        return 0
    else:
        logger.error(f"generate msdt fail")
        return -1
