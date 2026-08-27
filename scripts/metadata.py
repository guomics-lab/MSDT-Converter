import numpy as np
import pandas as pd
import os
import argparse
import shutil
import time
import json
import logging
import sqlite3
import xml.etree.ElementTree as ET

logger = logging.getLogger(__name__)

def d_metadata(bruker_d_folder_name):
    with sqlite3.connect(os.path.join(bruker_d_folder_name, "analysis.tdf")) as sql_database_connection:
        global_meta_data = pd.read_sql_query("SELECT * from GlobalMetaData",sql_database_connection)
        frames = pd.read_sql_query("SELECT * FROM Frames", sql_database_connection)
        if 9 in frames.MsMsType.values:
            acquisition_mode = "diaPASEF"
        elif 8 in frames.MsMsType.values:
            acquisition_mode = "ddaPASEF"
        else:
            acquisition_mode = "noPASEF"
    return(acquisition_mode, global_meta_data)

def get_d_metadata(bruker_d_folder_name):
    mode, metadata = d_metadata(bruker_d_folder_name)
    keys = metadata["Key"].tolist()
    metadata_list = metadata["Value"].tolist()
    meta_df = pd.DataFrame([metadata_list], columns=keys)
    meta_df["Filename"] = bruker_d_folder_name.split('/')[-1]
    meta_df["AcquisitionMode"] = mode
    meta_df["DatabaseSearchSoftware"] = "Sage"
    meta_df["DatabaseSearchSoftwareVersion"] = "0.14.7"
    REQUIRED_COLUMNS = ['Filename',
                       'AcquisitionMode',
                       'SchemaType',
                       'SchemaVersionMajor',
                       'SchemaVersionMinor',
                       'AcquisitionSoftwareVendor',
                       'InstrumentVendor',
                       'TimsCompressionType',
                       'ClosedProperly',
                       'MaxNumPeaksPerScan',
                       'AnalysisId',
                       'DigitizerNumSamples',
                       'PeakListIndexScaleFactor',
                       'MzAcqRangeLower',
                       'MzAcqRangeUpper',
                       'OneOverK0AcqRangeLower',
                       'OneOverK0AcqRangeUpper',
                       'AcquisitionSoftware',
                       'AcquisitionSoftwareVersion',
                       'AcquisitionFirmwareVersion',
                       'AcquisitionDateTime',
                       'InstrumentName',
                       'InstrumentFamily',
                       'InstrumentRevision',
                       'InstrumentSourceType',
                       'InstrumentSerialNumber',
                       'OperatorName',
                       'Description',
                       'SampleName',
                       'MethodName',
                       'DenoisingEnabled',
                       'DatabaseSearchSoftware',
                       'DatabaseSearchSoftwareVersion']
    missing_cols = [col for col in REQUIRED_COLUMNS if col not in meta_df.columns]
    for col in missing_cols:
        meta_df[col] = 'NotFound'
    meta_df = meta_df[REQUIRED_COLUMNS]
    return meta_df

def extract_filename(elem):
    file_name = "NotRawFile"
    name = elem.get('name')
    if name.endswith('.raw'):
        file_name = name
    return file_name

def extract_instrument(elem):
    instrument = "NotFound"
    user_param = elem.find('{http://psi.hupo.org/ms/mzml}userParam')
    if user_param is not None:
        value = user_param.get('value')
        if value != "":
            instrument = value
            return instrument
    else:
        for cv_param in elem.findall('{http://psi.hupo.org/ms/mzml}cvParam'):
            name = cv_param.get('name')
            value = cv_param.get('value')
            if value == "":
                instrument = name
                break
    return instrument

def extract_activation(elem):
    fragmentation_mode = "NotFound"
    collision_energy = "NotFound"
    for cv_param in elem.findall('{http://psi.hupo.org/ms/mzml}cvParam'):
        name = cv_param.get('name')
        value = cv_param.get('value')
        if value == "":
            fragmentation_mode = name
        if name == "collision energy":
            collision_energy = value
    try:
        collision_energy = float(collision_energy)
    except:
        collision_energy = "NotFound"
    return fragmentation_mode, str(collision_energy)

def extract_method(elem):
    method = "NotFound"
    user_param = elem.find('.//{http://psi.hupo.org/ms/mzml}userParam')
    name = user_param.get('name')
    if name != "":
        method = name
    return method

def get_mzml_metadata(file_path):
    context = ET.iterparse(file_path, events=('start', 'end'))
    _, root = next(context)
    result = {}
    found_tag = []
    all_tag = ['sourceFile','referenceableParamGroup','activation','dataProcessing']
    for event, elem in context:
        if event == 'end':
            tag_name = elem.tag.split('}')[-1]
            if tag_name == 'sourceFile' and tag_name not in found_tag:
                file_name = extract_filename(elem)
                result['Filename'] = file_name
                found_tag.append(tag_name)
            if tag_name == 'referenceableParamGroup' and tag_name not in found_tag:
                instrument = extract_instrument(elem)
                result['Instrument'] = instrument
                found_tag.append(tag_name)
            if tag_name == 'activation' and tag_name not in found_tag:
                fragmentation_mode, collision_energy = extract_activation(elem)
                result['FragmentationMode'] = fragmentation_mode
                result['CollisionEnergy'] = collision_energy
                found_tag.append(tag_name)
            if tag_name == 'dataProcessing' and tag_name not in found_tag:
                method = extract_method(elem)
                result['Method'] = method
                found_tag.append(tag_name)
            if set(found_tag) == set(all_tag):
                logger.info(f"metadata in {file_path} has been extracted")
                break
    root.clear()
    meta_df = pd.DataFrame([result])
    meta_df['DatabaseSearchSoftware'] = ['Sage']
    meta_df['DatabaseSearchSoftwareVersion'] = ['0.14.7']
    return meta_df

def extract_metadata(raw_path, save_path):
    """
    Extract metadata from raw mass spectrometry files.
    Supports Bruker .d (timsTOF) and .mzML formats.
    Return values:
        0: Successfully extracted
        -1: Extraction failed
    """
    try:
        if raw_path.endswith('.d'):
            meta_df = get_d_metadata(raw_path)
        elif raw_path.endswith('.mzML'):
            meta_df = get_mzml_metadata(raw_path)
        else:
            logger.error(f"unsupported file format: {raw_path}, only .d and .mzML are supported")
            return -1

        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        meta_df.to_csv(save_path, index=False)
        logger.info(f"metadata saved to {save_path}")
        return 0
    except Exception as e:
        logger.error(f"Error occurs when extracting metadata from {raw_path}: {e}")
        return -1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-raw', type=str, default="", help="path to raw spectrum file (.d or .mzML)")
    parser.add_argument('-save', type=str, default="", help="output path for metadata CSV")
    args = parser.parse_args()
    extract_metadata(args.raw, args.save)
