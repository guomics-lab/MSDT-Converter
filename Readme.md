# 📦 MSDT-Converter

**MassNet-Converter** is a tool for converting commonly used mass spectrometry data formats into the Mass Spectrometry DDA Tensor (MSDT) format—an efficient, standardized, and AI-friendly, 
schema-rich data representation designed for large-scale proteomics analysis, storage, exchange, and reuse.  

For more detailed information on the MSDT format, including its schema and design rationale, please refer to [MSDT v1.0_updated_20250211.pdf](https://github.com/guomics-lab/MSDT-Converter/blob/main/MSDT%20v1.0_updated_20250211.pdf) 

---

## ✨ Key Features

* **Supported Input Formats:**
    * `mzML` (standard open format)
    * `MGF` (Mascot Generic Format)
    * Bruker’s native `.d` directory format (TimsTOF)
    * WIFF-derived `mzML` (SCIEX, converted beforehand with a compatible tool)

* **Output Format:**
    * Standardized **MSDT** files stored in **Apache Parquet**, enabling efficient storage, high compression, and scalable data preparation for large-scale AI workflows.

* **Optimized for AI Workflows:**
    * Converts raw and search result data into structured **tensor format** for seamless integration with machine learning models, such as XuanjiNovo, DeepLC, and DDA-BERT.

* **Dockerized Deployment:**
    * Ready-to-use Docker image available on Docker Hub.
    * Run the converter in a reproducible environment without manual dependencies.

---

# 📥 Getting the Test Data and Configurations

We provide test datasets from **Thermo, SCIEX, and Bruker platforms, covering mzML, .d, and MGF formats**, together with configuration files for the FragPipe and Sage search engines and the corresponding converted MSDT files.

All test datasets, configuration files, and converted MSDT files can be downloaded from the Google Drive link below:

* **🔗 Download Link:** [test data and configs](https://drive.google.com/drive/folders/1gZnnue6BFuv6rowjXCX3NC_ug8cr2Pkn?usp=drive_link)

---

## 💻 Command Line Usage Examples (Docker)

Below are command line examples for running the data conversion using the locally built `msdt-converter:v2` Docker image for different instrument data. Build it first with `docker build -t msdt-converter:v2 .` from the repository root.

> **Note:** Please replace the local path `D:\Work\MSDT_Converter` in the commands with your actual data storage path.

### 1. Thermo Data Conversion

Uses the `config_mzml.json` configuration file.

```bash
docker run --rm -v "D:\Work\MSDT_Converter":/home msdt-converter:v2 -config=/home/config_mzml.json
```

### 2. Bruker Data Conversion

Uses the `config_tims.json` configuration file.

```bash
docker run --rm -v "D:\Work\MSDT_Converter":/home msdt-converter:v2 -config=/home/config_tims.json
```

### 3. SCIEX Data Conversion

Uses the `config_wiff.json` configuration file.

```bash
docker run --rm -v "D:\Work\MSDT_Converter":/home msdt-converter:v2 -config=/home/config_wiff.json
```


The running times below were measured on our test data using a single CPU core only, without GPU acceleration or multi-core parallelization. Actual running times may vary depending on input size, hardware and system load.

| | mzml | tims | wiff |
| :--- | :--- | :--- | :--- |
| **generate_rawspectrum** | 20 s | 3 min | 1 min |
| **generate_sage_search_result** | 30 s | 1 min | 1 min |
| **generate_fragpipe_search_result** | 7 min | 7 min | 5 min |
| **generate_msdt(sage)** | 5 s | 30 s | 1 min |
| **msdt_2_mgf** | 5 s | 1 min | 1 min |
| **convert_2_msdt** | 5 s | 5 s | 5 s |


---

# 🚀 Quick Start Guide

## Environment Setup

Setup is available via both Docker and Conda. Choose between Option A (Docker) and Option B (Conda) below.

## Option A: Docker

This repository provides a self-contained **Docker image** that encapsulates all necessary environments and
dependencies for the MassNet-DDA conversion utility. By using this image, users can quickly launch the tool without
complex setup.

### Prerequisites

* **Docker Desktop** (for Windows/Mac) or **Docker Engine** (for Linux) must be installed and running.

---

The process builds the image from this repository and then runs a container, mapping your local data
directory to the container's working directory. Docker automatically pulls the published v1.3 toolchain base image during the build.

1. **Ensure the Docker service is running.**
2. **Build the Docker Image** from the repository root:
   ```bash
   docker build -t msdt-converter:v2 .
   ```
3. **Run the Container** (Example using a typical Linux absolute path):
   ```bash
   docker run --rm -v /home/user/MassNet-DDA:/home/test_data msdt-converter:v2 -config=/home/test_data/config.json
   ```

## Option B: Conda
You can install MassNet-Converter in a Conda environment. This option is recommended if you prefer a Python-native setup or wish to modify the source code.

⏱️ Estimated setup time: **~2–5 minutes**

### Prerequisites

Download jdk11 from [here](https://guomics-share.oss-cn-shanghai.aliyuncs.com/SOFTWARE/MSDT-Converter/jdk-11.0.26.zip),
unzip and move to project root directory.
> **⚠️Note**: This project was built using **FragPipe v21.1**, which includes the following core components: **MSFragger v4.0**, **Philosopher v5.1.1**, **diaTracer v1.1.5**, **IonQuant v1.10.27**. If you're using a different version of FragPipe, make sure to download compatible versions of these tools and properly configure your environment to ensure smooth execution of the analysis pipeline.. 

You can download FragPipe from [here](https://guomics-share.oss-cn-shanghai.aliyuncs.com/SOFTWARE/MSDT-Converter/FragPipe-21.1.zip) with the following core components。


> **⚠️Note**: FragPipe versions may differ in their output fields, directory structures, and execution commands. When working with results generated by a different version, users are advised 
to consult the official output specifications for that release. Accordingly, configuration files (to be specified) or the parsing script (MSDT-Converter/scripts/search_engine.py) 
may need to be adjusted to ensure correct interpretation and processing of the data.

To learn more about FragPipe’s usage and configuration, please visit: https://github.com/Nesvilab/FragPipe.

1. Create a new conda environment first:

```
conda create --name msdt-converter python=3.13
```

This will create an anaconda environment

2. Activate this environment by running:

```
conda activate msdt-converter
```

3. Install dependencies:

```
pip install -r ./requirements.txt
```

4. Set up file permissions:

After cloning the repository and completing the installation, run the following command in the project’s root directory to ensure that all files and subdirectories have the appropriate access and execution permissions:

```
chmod -R 775 .
```

### Run the script
> **⚠️Note**: Before running the script, please download [test data and configs](https://drive.google.com/drive/folders/1gZnnue6BFuv6rowjXCX3NC_ug8cr2Pkn?usp=drive_link) to the root directory.

```bash
python convert.py enrich \
  --file-list "/test/run1/file.txt" \
  --data_type "mzml" \
  --workdir "/test/results" \
  --fasta "/test/Homo_sapiens_reviewed.fasta" \
  --workflow workflows/Default.workflow \
  --threads 10
```

---

## ⚙️ Configuration File (`config.json`)

The container requires a single **JSON configuration file** to define which steps to execute and to specify all
necessary input, output, and processing parameters.

---

### 📑 Overall Workflow Steps

The configuration is structured by the main processing steps. Each primary object controls a specific function.

| Parameter Name | Description |
| :--- | :--- |
| **`generate_rawspectrum`** | Parameters for extracting raw spectral data into a `.tsv` file. |
| **`generate_sage_search_result`** | Parameters for running the **Sage** search engine. |
| **`generate_fragpipe_search_result`** | Parameters for running the **FragPipe** search pipeline. |
| **`generate_msdt`** | Parameters for converting search results and raw data into the **MSDT** format. |
| **`convert_2_msdt`** | Parameters for converting other formats (like MGF) directly to MSDT. |
| **`msdt_2_mgf`** | Parameters for converting MSDT back to the MGF format. |

---

### 1️⃣ `generate_rawspectrum`

| Parameter | Data Type | Example Value | Description |
| :--- | :--- | :--- | :--- |
| **`need`** | `boolean` | `true` | Set to `true` to execute this step (extract raw spectra). |
| **`data_type`** | `string` | `"mzml"` | The type of input data: `mzml`, `tims`, or `wiff2mzml` (for mzML converted from WIFF). |
| **`data_path`** | `string` | `/home/test_data/.../DDA_ingel_3D.mzML` | **Input.** Absolute path to the raw data file (relative to the Docker mounted volume). |
| **`output`** | `string` | `/home/test_data/.../3D_rawspectrum.tsv` | **Output.** Path for the generated raw spectrum TSV file. |

---

### 2️⃣ `generate_sage_search_result`

| Parameter | Data Type | Example Value | Description |
| :--- | :--- | :--- | :--- |
| **`need`** | `boolean` | `true` | Set to `true` to execute this step (run Sage search). |
| **`workdir`** | `string` | `/home/test_data/2_generate_sage_search_result` | **Input/Output.** Working directory where Sage will generate its result files. |
| **`fasta`** | `string` | `/home/test_data/.../Homo_sapiens_reviewed.fasta` | **Input.** Path to the FASTA protein sequence database file. |
| **`data_path`** | `string` | `/home/test_data/.../DDA_ingel_3D.mzML` | **Input.** Path to the mzML file used for searching. |
| **`config_path`** | `string` | `/home/test_data/.../sage_config.json` | **Input.** Path to the specific configuration file for the Sage search engine. |

---

### 3️⃣ `generate_fragpipe_search_result`

| Parameter | Data Type | Example Value | Description |
| :--- | :--- | :--- | :--- |
| **`need`** | `boolean` | `true` | Set to `true` to execute this step (run FragPipe search). |
| **`workdir`** | `string` | `/home/test_data/3_generate_fragpipe_search_result` | **Input/Output.** Working directory where FragPipe will generate results. |
| **`data_path`** | `string` | `/home/test_data/.../DDA_ingel_3D.mzML` | **Input.** Path to the mzML file used for searching. |
| **`workflow_path`** | `string` | `/home/test_data/.../LFQ_DDA_human_noNQ.workflow` | **Input.** Path to the FragPipe workflow configuration file and the fasta path should be set in the workflow. |
| **`manifest_path`** | `string` | `/home/test_data/.../fragpipe-files.fp-manifest` | **Output.** Path for the FragPipe temporary manifest output file. |
| **`thread_num`** | `integer` | `10` | The number of CPU threads to use for the FragPipe search process. |

---

### 4️⃣ `generate_msdt`

This section contains nested configurations based on data type (`tims`, `mzml`, `wiff`).

#### **4.1. `generate_msdt` -> `tims`**

| Parameter | Data Type | Example Value | Description |
| :--- | :--- | :--- | :--- |
| **`need_tims`** | `boolean` | `false` | Set to `true` to generate MSDT from **tims** data (not currently configured in the example). |
| **`rawspectrum_path`** | `string` | `""` | **Input.** Path to the raw spectrum file. |
| **`sage_search_result_path`** | `string` | `""` | **Input.** Path to the Sage search result file. |
| **`unify_residue`** | `boolean` | `true` | If `true`, the residue format will be converted to the unified MSDT format. |
| **`output`** | `string` | `""` | **Output.** Path for the generated Sage MSDT file. |

#### **4.2. `generate_msdt` -> `mzml`**

| Parameter | Data Type | Example Value | Description |
| :--- | :--- | :--- | :--- |
| **`need_mzml`** | `boolean` | `true` | Set to `true` to generate MSDT from **mzML** related data. |
| **`need_sage`** | `boolean` | `true` | Set to `true` to generate MSDT from **Sage** search results. |
| **`need_fragpipe`** | `boolean` | `true` | Set to `true` to generate MSDT from **FragPipe** search results. |
| **`rawspectrum_path`** | `string` | `/home/test_data/.../3D_rawspectrum.tsv` | **Input.** Path to the raw spectrum file. |
| **`sage_search_result_path`** | `string` | `/home/test_data/.../D_search_result.tsv` | **Input.** Path to the Sage search result file. |
| **`fp_pin_path`** | `string` | `/home/test_data/.../A18..._edited.pin` | **Input.** Path to the FragPipe `.pin` file. |
| **`percolator_target_path`** | `string` | `/home/test_data/.../exp/..._target_psms.tsv` | **Input.** Percolator target PSM TSV from the same FragPipe result directory as the PIN. |
| **`percolator_decoy_path`** | `string` | `/home/test_data/.../exp/..._decoy_psms.tsv` | **Input.** Percolator decoy PSM TSV from the same FragPipe result directory as the PIN. |
| **`run_id`** | `string` | `"sample01"` | Use only for pooled global-Percolator TSVs whose PSM IDs start with `sample01::`; omit for ordinary per-run TSVs. |
| **`fdr_threshold`** | `number` | `0.01` | Optional q-value threshold for targets; all decoys are retained. |
| **`sage_unify_residue`** | `boolean` | `true` | If `true`, Sage residue format converts to MSDT format. |
| **`fp_unify_residue`** | `boolean` | `true` | If `true`, FragPipe residue format converts to MSDT format. |
| **`sage_output`** | `string` | `/home/test_data/.../sage_msdt.parquet` | **Output.** Path for the generated Sage MSDT `.parquet` file. |
| **`fp_output`** | `string` | `/home/test_data/.../fp_msdt.parquet` | **Output.** Path for the generated FragPipe MSDT `.parquet` file. |

#### **4.3. `generate_msdt` -> `wiff`**

| Parameter | Data Type | Example Value | Description |
| :--- | :--- | :--- | :--- |
| **`need_wiff`** | `boolean` | `false` | Set to `true` to generate MSDT from **WIFF** related data (not currently configured in the example). |
| **`need_sage`** | `boolean` | `true` | Set to `true` to generate a Sage-derived WIFF MSDT. |
| **`need_fragpipe`** | `boolean` | `true` | Set to `true` to generate a FragPipe/Percolator-derived WIFF MSDT. |
| **`wiff_mzml_path`** | `string` | `""` | **Input.** Path to the mzML file converted from WIFF. |
| **`rawspectrum_path`** | `string` | `""` | **Input.** Path to the raw spectrum file. |
| **`sage_search_result_path`** | `string` | `""` | **Input.** Path to the Sage search result file. |
| **`fp_pin_path`** | `string` | `""` | **Input.** FragPipe edited PIN used to build the FP MSDT rows. |
| **`percolator_target_path`** | `string` | `""` | **Input.** Percolator target PSM TSV retained by the workflow. |
| **`percolator_decoy_path`** | `string` | `""` | **Input.** Percolator decoy PSM TSV retained by the workflow. |
| **`run_id`** | `string` | `"sample01"` | Use only for pooled global-Percolator TSVs whose PSM IDs start with `sample01::`; omit for ordinary per-run TSVs. |
| **`fdr_threshold`** | `number` | `0.01` | Optional q-value threshold for targets; all decoys are retained. |
| **`sage_unify_residue`** | `boolean` | `true` | Normalize Sage modified sequences to the MSDT representation. |
| **`fp_unify_residue`** | `boolean` | `true` | Normalize FragPipe/Percolator modified sequences before PSM matching. |
| **`sage_output`** | `string` | `""` | **Output.** Path for the generated Sage MSDT file. |
| **`fp_output`** | `string` | `""` | **Output.** Path for the FP-derived MSDT containing `score`, `q-value`, and `PEP`. |

---

### 5️⃣ `convert_2_msdt`

This section handles direct conversion from other data formats to MSDT.

#### **5.1. `convert_2_msdt` -> `mgf`**

| Parameter | Data Type | Example Value | Description |
| :--- | :--- | :--- | :--- |
| **`need`** | `boolean` | `true` | Set to `true` to execute this MGF conversion step. |
| **`mgf_path`** | `string` | `/home/test_data/.../180624_G12.MGF` | **Input.** Path to the MGF file to be converted. |
| **`output_path`** | `string` | `/home/test_data/.../180624_G12.parquet` | **Output.** Path for the generated MSDT `.parquet` file. |
| **`field_type_dict`** | `object` | `{...}` | A dictionary defining the fields present in the MGF file and their corresponding data types. |

> **`field_type_dict`** details:

| Key | Data Type | Description |
| :--- | :--- | :--- |
| `TITLE` | `"string"` | The title of the spectrum (required). |
| `PEPMASS` | `"float"` | The precursor mass (required). |
| `CHARGE` | `"int"` | The precursor charge (e.g., "2+"). Must be convertible to integer. |
| `RTINSECONDS` | `"float"` | The retention time in seconds. |
| `INSTRUMENT` | `"string"` | The instrument name. |

---

### 6️⃣ `msdt_2_mgf`

| Parameter | Data Type | Example Value | Description |
| :--- | :--- | :--- | :--- |
| **`need`** | `boolean` | `true` | Set to `true` to execute this step (convert MSDT back to MGF). |
| **`msdt_path`** | `string` | `/home/test_data/.../sage_msdt.parquet` | **Input.** Path to the MSDT `.parquet` file to be converted. |
| **`output_path`** | `string` | `/home/test_data/.../sage.mgf` | **Output.** Path for the generated MGF file. |

---

## 📋 Extract Raw File Metadata (`metadata.py`)

The `scripts/metadata.py` script extracts acquisition metadata from raw mass spectrometry files and exports it as a CSV file.  
**Supported formats:** `.mzML` (standard open format) and `.d` (Bruker TimsTOF native directory).

### Usage

```bash
python scripts/metadata.py -raw <input.d or input.mzML> -save <output.csv>
```

| Parameter | Description |
| :--- | :--- |
| **`-raw`** | Path to the raw file: a Bruker `.d` directory or an `.mzML` file. |
| **`-save`** | Output path for the metadata CSV file. |

### Example

```bash
# Extract metadata from a Bruker TimsTOF .d folder
python scripts/metadata.py -raw /data/sample.d -save /data/sample_metadata.csv

# Extract metadata from an mzML file
python scripts/metadata.py -raw /data/sample.mzML -save /data/sample_metadata.csv
```

### Output Fields

For Bruker `.d` files, the output includes 33 standardized metadata fields such as `Filename`, `AcquisitionMode` (diaPASEF / ddaPASEF / noPASEF), `InstrumentName`, `AcquisitionDateTime`, `SampleName`, `MethodName`, and more.

For `.mzML` files, the output includes `Filename`, `Instrument`, `FragmentationMode`, `CollisionEnergy`, `Method`, and database search software version info.

---

## 📚 Citation

When using **MassNet-Converter**, please cite the following publication:

Jun, A., Zhang, X., Zhang, X., Wei, J., Zhang, T., Deng, Y., ... & Guo, T. (2025). MassNet: billion-scale AI-friendly mass spectral corpus enables robust de novo peptide sequencing. bioRxiv, 2025-06.
