# PRIDE Real-Data Test Candidates for MSDT-Converter v2

> Research date: 2026-08-11. This review used only official PRIDE Archive and ProteomeXchange project metadata, the PRIDE API, and official file hosts; no large files were downloaded. File sizes are based on the API `fileSizeBytes` value. MiB = bytes / 1,048,576.

## Recommendation

The best choice is a three-tier test set rather than one large project:

1. **Primary SCIEX smoke and batch test: `PXD061973`.** One `WIFF/WIFF.scan` pair is only 144.81 MiB, while all four pairs total 615.28 MiB. This dataset can validate WIFF scan mapping, `file_list`, and global FDR together.
2. **Thermo regression test: `PXD000001`.** The classic LTQ Orbitrap Velos RAW file is only 210.26 MiB. The project also provides a FASTA file, a Mascot DAT file, and a later mzML file linked directly from the project record, making it suitable for RAW-to-mzML conversion regression testing.
3. **PTM and multiple-charge stress test: `PXD079474`.** The converted mzML, FASTA, and mzIdentML files total approximately 587 MiB. The official method specifies Exploris 480 DDA Top15 acquisition, HCD, precursor charges 2-6, and an Fe-NTA-enriched STY phosphopeptide sample.

For an official one-to-one reference consisting of an mzML file and a same-name WIFF file, add `PXD064530`. It is valuable, but a complete run requires 2.44 GiB and is therefore not recommended for the first smoke test.

## A. Preferred: PXD061973 (Small SCIEX Files and Four-File Batch)

PRIDE project: [PXD061973](https://www.ebi.ac.uk/pride/archive/projects/PXD061973)

Official metadata: [project API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD061973); [files API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD061973/files?pageSize=1000&page=0)

- Instrument: the PRIDE controlled field specifies **TripleTOF 5600**, while the sample-method text specifies **Triple-TOF 5600+ (AB SCIEX)**.
- Acquisition: the PRIDE controlled field specifies **data-dependent acquisition**, while the method text specifies **information-dependent acquisition (IDA)**.
- PTMs: PRIDE states that no PTMs are included in the dataset, which makes the project suitable for an initial test without PTM-related complexity.
- Processing: the original authors used ProteinPilot 4.5, an E. coli plus rabbit UniProt database dated 2021-09-28, and 1% FDR. However, this is a partial submission and **does not include a FASTA file or search results**. It is therefore suitable for a standardized FragPipe re-search with this project, but not for reproducing the identification count in the original paper.

| Run | WIFF | WIFF.scan | Total |
|---|---:|---:|---:|
| `1` | `1.wiff` — 12,455,936 B (11.88 MiB) | `1.wiff.scan` — 139,391,092 B (132.93 MiB) | 144.81 MiB |
| `2` | `2.wiff` — 12,156,928 B (11.59 MiB) | `2.wiff.scan` — 155,982,256 B (148.76 MiB) | 160.35 MiB |
| `3` | `3.wiff` — 12,161,024 B (11.60 MiB) | `3.wiff.scan` — 149,001,820 B (142.10 MiB) | 153.70 MiB |
| `4` | `4.wiff` — 12,840,960 B (12.25 MiB) | `4.wiff.scan` — 151,177,880 B (144.17 MiB) | 156.42 MiB |
| **All runs** |  |  | **645,167,896 B = 615.28 MiB** |

Official HTTPS directory prefix: `https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/03/PXD061973/`. Examples: [1.wiff](https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/03/PXD061973/1.wiff); [1.wiff.scan](https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/03/PXD061973/1.wiff.scan).

**Usage:** Download only Run 1 for the initial smoke test. After it passes, download Runs 1-4, convert all four to mzML with identical MSConvert parameters, search them with the same FASTA and workflow, and then test `file_list` and global Percolator. **A `.wiff` file and its same-name `.wiff.scan` companion must be downloaded and retained together.**

## B. Thermo Regression: PXD000001

PRIDE project: [PXD000001](https://www.ebi.ac.uk/pride/archive/projects/PXD000001)

Official metadata: [project API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD000001); [files API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD000001/files?pageSize=1000&page=0)

- Instrument: **LTQ Orbitrap Velos**. The experiment type is bottom-up proteomics. The original filename contains `Top10HCD`, making it useful for HCD regression testing. The PRIDE controlled field does not separately specify DDA, so the acquisition mode should be confirmed from the actual mzML scan structure.
- PTMs: PRIDE lists monohydroxylation, TMT6plex acylation, and methylthiolation.

| Role | File | Size |
|---|---|---:|
| Thermo RAW | [`TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw) | 220,475,548 B (210.26 MiB) |
| Database | [`erwinia_carotovora.fasta`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/erwinia_carotovora.fasta) | 1,657,668 B (1.58 MiB) |
| Original search result (reference only) | [`F063721.dat`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/F063721.dat) | 21,185,462 B (20.20 MiB) |
| Later converted mzML | [`...01-20141210.mzML`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML) | 450,032,788 B (429.18 MiB) |

Note: The current files API lists the RAW, FASTA, DAT, and an older mzXML file, but does not list the `-20141210.mzML` file shown above. However, the project's own data-processing protocol links directly to this mzML file, and the official PRIDE HTTPS host returns status 200 with the Content-Length shown above. The file can therefore serve as a conversion reference, but it should not be assumed to use exactly the same MSConvert parameters as the current workflow.

**Minimum download:** RAW + FASTA = 211.84 MiB.

**RAW/mzML comparison:** RAW + FASTA + official mzML = 641.03 MiB.

## C. PTM and Multiple-Charge Stress Test: PXD079474

PRIDE project: [PXD079474](https://www.ebi.ac.uk/pride/archive/projects/PXD079474)

Official metadata: [project API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD079474); [files API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD079474/files?pageSize=1000&page=0)

- Instrument: **Orbitrap Exploris 480**.
- Sample: Fe-NTA-enriched phosphopeptides. Variable modifications in the search include oxidation (M), protein N-terminal acetylation, and phosphorylation (S/T/Y), with fixed carbamidomethylation (C). The PRIDE controlled PTM field is phosphorylated residue.
- Acquisition: the method text explicitly specifies **DDA Top15, HCD, and MS2 charge 2-6**. This dataset can directly stress-test `charge`, modified sequences, and PSM enrichment fields. An allowed charge range of 2-6 does not guarantee that every charge occurs; the actual distribution must be verified after extraction.

| Role | File | Size |
|---|---|---:|
| Preferred input | [`SpermC_Rep1.mzML`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD079474/SpermC_Rep1.mzML) | 591,044,161 B (563.66 MiB) |
| FASTA | [`SwissProt_Mouse_03232022.fasta`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD079474/SwissProt_Mouse_03232022.fasta) | 11,926,463 B (11.37 MiB) |
| External identification reference | [`SpermC_Rep1.mzid`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD079474/SpermC_Rep1.mzid) | 12,712,011 B (12.12 MiB) |
| Optional Thermo RAW | [`SpermC_Rep1.raw`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD079474/SpermC_Rep1.raw) | 1,268,241,220 B (1,209.49 MiB) |
| Optional PD result | [`SpermC_Rep1.pdResult`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD079474/SpermC_Rep1.pdResult) | 593,108,992 B (565.63 MiB) |

**Recommended download:** Start with the mzML, FASTA, and mzIdentML files, totaling 587.16 MiB. This is sufficient for PTM and charge checks. There is no need to download the 1.21 GiB RAW file or the 565.63 MiB `pdResult` file unless Exploris RAW conversion will also be tested.

## D. Optional Gold Standard: PXD064530 (Same-Name WIFF/WIFF.scan/mzML)

PRIDE project: [PXD064530](https://www.ebi.ac.uk/pride/archive/projects/PXD064530)

Official metadata: [project API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD064530); [files API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD064530/files?pageSize=1000&page=0)

- Instrument and acquisition: **TripleTOF 5600, DDA**.
- PTMs: PRIDE lists oxidation, deamidation, iodoacetamide derivatization, and iTRAQ8plex-116 acylation. The processing method describes search settings for iTRAQ8plex (K/Y/peptide N-terminus), carbamidomethyl (C), and other modifications in detail.
- Primary value: the project's data-processing protocol explicitly states that the `.wiff` file was converted to mzML using the ProteoWizard `qtofpeakpicker` with resolution 15,000 and threshold 7.5, and all three files share the same base name. This makes the project particularly suitable for validating native IDs, MS2 counts, retention time, and precursor-m/z mapping.

| File | Size |
|---|---:|
| [`23062_EA1Cleavage_iTRAQ.wiff`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/11/PXD064530/23062_EA1Cleavage_iTRAQ.wiff) | 10,194,944 B (9.72 MiB) |
| [`23062_EA1Cleavage_iTRAQ.wiff.scan`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/11/PXD064530/23062_EA1Cleavage_iTRAQ.wiff.scan) | 1,551,030,628 B (1,479.18 MiB) |
| [`23062_EA1Cleavage_iTRAQ.mzML`](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/11/PXD064530/23062_EA1Cleavage_iTRAQ.mzML) | 1,062,876,923 B (1,013.64 MiB) |
| **Total** | **2,624,102,495 B = 2.44 GiB** |

## E. Optional SCIEX and MGF Reference: PXD074970

PRIDE project: [PXD074970](https://www.ebi.ac.uk/pride/archive/projects/PXD074970)

Official metadata: [project API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD074970); [files API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD074970/files?pageSize=1000&page=0)

Each run in this project provides `WIFF + WIFF.scan + MGF`, so the MGF file can serve as an independent peak-list reference. The instrument is a TripleTOF 5600, and the ProteinPilot/Paragon search method includes variable Met oxidation, N-terminal pyroGlu/acetylation, and Cys carboxyamidomethylation.

| Run | WIFF | WIFF.scan | MGF | Total |
|---|---:|---:|---:|---:|
| `pcDNA_section_1` | 17.94 MiB | 474.21 MiB | 2.85 MiB | 495.00 MiB |
| `PACS1_WT_section_2` | 16.78 MiB | 461.44 MiB | 7.31 MiB | 485.53 MiB |
| `PACS1_R203W_section_3` | 16.62 MiB | 449.72 MiB | 5.33 MiB | 471.67 MiB |

Exact filenames are available from the [files API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD074970/files?pageSize=1000&page=0). One metadata detail requires caution: the controlled experiment-type field says “Top-down proteomics,” while the sample and data-processing text explicitly describes trypsin and peptide-based ProteinPilot identification. The project is suitable for format and peak-list comparisons, but a test report should not repeat the “Top-down” label without qualification.

## Recommended Execution Order and Acceptance Criteria

### 1. Minimum SCIEX Smoke Test (~145 MiB)

Use only `PXD061973/1.wiff` and `1.wiff.scan`:

- Convert to mzML on Windows with a SCIEX-compatible MSConvert environment.
- Confirm that `wiff_mzml_rawspecturm` extracts all MS2 spectra.
- Confirm that FragPipe `ScanNr` maps to Parquet `scan + 1` without missing or duplicate matches.
- Do not use the non-unique cycle value alone to align native IDs.
- Require all retention-time and precursor-m/z cross-checks to pass.

### 2. Four-Run `file_list` and True Global FDR (~615 MiB)

Use Runs 1-4 from `PXD061973`, with exactly the same FASTA, FragPipe workflow, and Percolator version for every run:

- Confirm that all four `file_list` inputs are actually searched and that outputs are not overwritten.
- Confirm that the feature headers and `DefaultDirection` values are consistent across all four PIN files.
- After global scan remapping, require `(run_id, ScanNr)` to be globally unique while preserving a shared remapped scan for multiple candidates from the same spectrum.
- After mapping global target and decoy results back to individual runs, confirm that PSM counts agree with the expected FDR filtering and that cross-run scan-number collisions do not discard data.

### 3. Thermo RAW Regression (Minimum ~212 MiB)

Use the `PXD000001` RAW and FASTA files. Compare schema, scan, retention time, precursor m/z, charge, and `score/q-value/PEP` missingness against the existing Thermo test results. If the later official mzML file is also downloaded, compare the MS2 count and key metadata. Explainable differences caused by conversion-tool or parameter changes are acceptable.

### 4. Phosphorylation and Multiple Charges (~587 MiB)

Use the mzML, FASTA, and mzIdentML files from `PXD079474`:

- Calculate the observed `charge` distribution and confirm that it is compatible with the DDA charge range of 2-6.
- Confirm that phosphorylation on S/T/Y, oxidation on M, N-terminal acetylation, and carbamidomethylation on C are represented without loss, duplication, or incorrect normalization.
- Compare FragPipe/Percolator output with the archived mzIdentML for scan, charge, and sequence coverage, but do not require one-to-one PSM agreement between different search engines.

## Input Types Not Recommended for the First Test Round

- A `.wiff` file without its `.wiff.scan` companion.
- Multi-gigabyte SWATH/DIA queues at the beginning of testing, because the current acceptance target is the DDA plus FragPipe/Percolator path.
- PIN files from different batches, FASTA databases, or search spaces forced into one global FDR analysis.
- Archived Mascot DAT, ProteinPilot XLSX, or PD `pdResult` files treated as FragPipe PIN input for this converter. These files are external references only.

## Official Resources

- [PRIDE Archive API guide](https://www.ebi.ac.uk/pride/ws/archive/v2/docs/api-guide.html): endpoints for projects, files, MSRun metadata, and related records.
- [PRIDE API overview](https://www.ebi.ac.uk/pride/markdownpage/prideapi): v3 project and search endpoints plus the Swagger entry point.
- [PRIDE file download guide](https://www.ebi.ac.uk/pride/markdownpage/pridefiledownload): FTP, Aspera, Globus, and streaming download methods.
- [ProteomeXchange dataset lookup](https://proteomecentral.proteomexchange.org/cgi/GetDataset): look up a ProteomeXchange record by PXD accession.
