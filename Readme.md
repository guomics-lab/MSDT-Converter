# 📦 MSDT-Converter

**MassNet-Converter** is a tool for converting commonly used mass spectrometry data formats into the Mass Spectrometry DDA Tensor (MSDT) 
format—an efficient, standardized, and AI-friendly representation designed for large-scale proteomics analysis.


✨ **Key Features**

**Supported input formats:**

* **mzML (standard open format)**

* **MGF (Mascot Generic Format)**

* **Bruker’s native .d format (TimsTOF)**

**Output format:**

Standardized MSDT files stored in **Apache Parquet**, enabling fast I/O, high compression, and compatibility with distributed GPU/TPU training pipelines.

**Optimized for AI workflows:**
Converts raw and search result data into structured tensor format for seamless integration with machine learning models, such as DeepLC, XuanjiNovo, and DDA-BERT.

**Dockerized deployment:**

* **Ready-to-use Docker image available on Docker Hub**

* **Run the converter in a reproducible environment without manual dependencies**

📄 **Citation**

If you use MassNet-Converter in your work, please cite:

Jun, A., Zhang, X., Zhang, X., Wei, J., Zhang, T., Deng, Y., ... & Guo, T. (2025). MassNet: billion-scale AI-friendly mass spectral corpus enables robust de novo peptide sequencing. bioRxiv, 2025-06.


# 🚀 Quick Start Guide

## Environment Setup

We provide both Docker and Conda set-up guide, user can choose between option A: Docker and option B: Conda below:

## Option A: Docker

### Prerequisites

* **Docker Desktop** (for Windows/Mac) or **Docker Engine** (for Linux) must be installed and running.

---

### 💻 1. Windows

The process involves **pulling the image from Docker Hub** and then running a container, mapping your local data
directory to the container's working directory.

1. **Open Docker Desktop.** Ensure the Docker engine is running.
2. **Pull the Docker Image** from the registry using your command line (e.g., PowerShell or Command Prompt):
   ```bash
   docker pull guomics2017/msdt-converter:v1.1
   ```
3. **Run the Container** by mounting your local working directory (`D:\Work\MassNet-DDA` in this example) to the
   container's internal data path (`/home/test_data`) and specifying the path to your configuration file (
   `config.json`):
   ```bash
   docker run --rm -v "D:\Work\MassNet-DDA":/home/test_data guomics2017/msdt-converter:v1.1 -config=/home/test_data/config.json
   ```
    * **Note:** The `-v` flag maps your local directory to the container. The paths must be adjusted according to your
      actual setup.

---

### 🐧 2. Linux

The process involves **pulling the image from Docker Hub** and then running a container, mapping your local data
directory to the container's working directory.

1. **Ensure the Docker service is running.**
2. **Pull the Docker Image** from the registry in your terminal:
   ```bash
   docker pull guomics2017/msdt-converter:v1.1
   ```
3. **Run the Container** (Example using a typical Linux absolute path):
   ```bash
   docker run --rm -v /home/user/MassNet-DDA:/home/test_data guomics2017/msdt-converter:v1.1 -config=/home/test_data/config.json
   ```

---

### 🍎 3. macOS

The process involves **pulling the image from Docker Hub** and then running a container, mapping your local data
directory to the container's working directory.

1. **Open Docker Desktop.** Ensure the Docker engine is running.
2. **Pull the Docker Image** from the registry in your terminal:
   ```bash
   docker pull guomics2017/msdt-converter:v1.1
   ```
3. **Run the Container** (Example using a typical macOS path):
   ```bash
   docker run --rm -v /Users/yourname/Documents/MassNet-DDA:/home/test_data guomics2017/msdt-converter:v1.1 -config=/home/test_data/config.json
   ```

## Option B: Conda

### Prerequisites

Download jdk11 from [here](https://guomics-share.oss-cn-shanghai.aliyuncs.com/SOFTWARE/MSDT-Converter/jdk-11.0.26.zip),
unzip and move to project root directory.
> **⚠️Note**: The jdk is from Oracle, we only provide the jdk for download easy.

Download FragPipe
from [here](https://guomics-share.oss-cn-shanghai.aliyuncs.com/SOFTWARE/MSDT-Converter/FragPipe-21.1.zip), unzip and
move to project root directory.
> **⚠️Note**: The FragPipe is from https://github.com/Nesvilab/FragPipe, we only provide the FragPipe for download easy.
The version contain some plugins.

Create a new conda environment first:

```
conda create --name msdt-converter python=3.13
```

This will create an anaconda environment

Activate this environment by running:

```
conda activate msdt-converter
```

then install dependencies:

```
pip install -r ./requirements.txt
```

### Run the script:

```bash
python convert.py -config=/home/test_data/config.json
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
| **`fasta_path`** | `string` | `/home/test_data/.../fasta.fas` | **Input.** Path to the FASTA protein sequence database file. |
| **`workflow_path`** | `string` | `/home/test_data/.../LFQ_DDA_human_noNQ.workflow` | **Input.** Path to the FragPipe workflow configuration file. |
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
| **`sage_unify_residue`** | `boolean` | `true` | If `true`, Sage residue format converts to MSDT format. |
| **`fp_unify_residue`** | `boolean` | `true` | If `true`, FragPipe residue format converts to MSDT format. |
| **`sage_output`** | `string` | `/home/test_data/.../sage_msdt.parquet` | **Output.** Path for the generated Sage MSDT `.parquet` file. |
| **`fp_output`** | `string` | `/home/test_data/.../fp_msdt.parquet` | **Output.** Path for the generated FragPipe MSDT `.parquet` file. |

#### **4.3. `generate_msdt` -> `wiff`**

| Parameter | Data Type | Example Value | Description |
| :--- | :--- | :--- | :--- |
| **`need_wiff`** | `boolean` | `false` | Set to `true` to generate MSDT from **WIFF** related data (not currently configured in the example). |
| **`wiff_mzml_path`** | `string` | `""` | **Input.** Path to the mzML file converted from WIFF. |
| **`rawspectrum_path`** | `string` | `""` | **Input.** Path to the raw spectrum file. |
| **`sage_search_result_path`** | `string` | `""` | **Input.** Path to the Sage search result file. |
| **`unify_residue`** | `boolean` | `true` | If `true`, the residue format will be converted to the unified MSDT format. |
| **`output`** | `string` | `""` | **Output.** Path for the generated Sage MSDT file. |

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