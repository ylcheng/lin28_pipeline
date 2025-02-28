# Lin28 Lineage Barcode Analysis Pipeline

## Overview
The **Lin28 Pipeline** processes **lineage barcode sequencing data** to identify cell lineages, correct sequencing errors, and analyze clonal relationships. It enables filtering, clustering, and correction of lineage barcodes using **Hamming distance-based filtering** and **Starcode clustering**.

## Features
- **Extract lineage barcodes** from FASTQ files
- **Filter and correct errors** using Hamming distance constraints
- **Cluster barcodes** using Starcode for lineage correction
- **Analyze clone relationships** based on barcode similarity
- **Generate outputs** for downstream lineage analysis

## Requirements

- **Python 3.x**
- Required Python libraries: `argparse`, `subprocess`, `logging`, `json`, `datetime`, `os`, `re`, `bz2`, `pandas`, `numpy`
- **Starcode** (for lineage barcode clustering)

## Usage
### **Basic Command**
```sh
python /path/to/lin28_pipeline.py \
    --code_path /path/to/lin28_pipeline \
    --folder_path /path/to/output_directory \
    --prefix sample_prefix \
    --fastq_pattern "*001.fastq.gz" \
    --max_hamming_distance 1 \
    --max_lineage_mismatch 0 \
    --max_start_flanking_mismatch 1 \
    --max_stop_flanking_mismatch 1 \
    --rna_file /path/to/rna_dataset.txt \
    --segment_match \
    --generate_clones_associated_files_before_correction \
    --compute_hamming_distances_within_and_between_clones \
    --starcode_path /path/to/starcode \
    --cluster_method sphere \
    --levenshtein_distance 2 \
    --max_correction 2
```

### **Command-line Arguments**
| Argument | Description |
|----------|-------------|
| `--code_path` | Path to the directory containing the pipeline scripts (default: current directory) |
| `--folder_path` | Input/output directory for processing (required) |
| `--prefix` | Base prefix for output file naming (required) |
| `--fastq_pattern` | Pattern to search for FASTQ files (default: `"*fastq.gz"`) |
| `--max_hamming_distance` | Maximum Hamming distance allowed for cell barcode matching (default: `0`) |
| `--max_lineage_mismatch` | Maximum mismatches allowed for lineage barcodes (default: `0`) |
| `--max_start_flanking_mismatch` | Maximum allowable mismatches in the start flanking sequence (default: `1`) |
| `--max_stop_flanking_mismatch` | Maximum allowable mismatches in the stop flanking sequence (default: `1`) |
| `--rna_file` | Path to the RNA dataset file (required). The file must be a tab-separated (`.tsv` or `.txt`) text file with **no header**. The first column must contain cell barcodes, and additional columns (if any) will be ignored. |
| `--segment_match` | Enable segment match filtering for cell barcode matching (default: `False`) |
| `--generate_clones_associated_files_before_correction` | Generate lineage per cell barcode and clone information files before barcode correction (default: `False`) |
| `--compute_hamming_distances_within_and_between_clones` | Compute Hamming distances within and between clones (default: `False`) |
| `--starcode_path` | Path to the Starcode executable (required) |
| `--cluster_method` | Clustering method (`sphere` or `message_passing`) (default: `sphere`) |
| `--levenshtein_distance` | Levenshtein distance threshold for clustering (default: `2`) |
| `--max_correction` | Maximum allowable nucleotide corrections for lineage barcodes (default: `2`) |

## Pipeline Workflow
1. **Extract Reads with Flanking Sequences** – Identify reads with known start/stop sequences.
2. **Extract Lineage Barcodes** – Parse FASTQ files to extract barcode information.
3. **Compute Hamming Distance** – Compare cell barcodes between RNA and lineage datasets.
4. **Filter & Correct Barcodes** – Apply mismatch constraints to filter erroneous reads.
5. **Collapse Reads Per UMI** – Count reads per UMI.
6. **Collapse UMIs Per Barcode** – Group lineage barcodes per cell barcode.
7. **Compute Clone Similarities** – Measure barcode similarity within and between clones.
8. **Starcode Clustering** – Perform lineage barcode correction with Starcode.
9. **Assign Clone IDs** – Cluster cells into clones based on barcode similarity.

## Output Files
| **File Name Pattern** | **Description** |
|----------------------|----------------|
| `*.lineage_barcode_reads.tsv` | Extracted lineage barcode reads from FASTQ files |
| `*.ambiguous_lineage_barcode_reads.tsv` | Reads with uncertain lineage barcode assignments |
| `*.filtered_lineage_barcode_reads.tsv` | Lineage barcode reads after filtering mismatches |
| `*.reads_per_umi.tsv` | Reads collapsed per UMI |
| `*.umi_per_lineage.tsv` | UMIs collapsed per lineage barcode |
| `*.lineages_per_cellbc.tsv` | Summary of lineage barcodes per cell barcode |
| `*.clone_info.tsv` | Clone grouping information based on lineage similarity |
| `*.lineage_hammings_within_and_btw_clones.tsv.bz2` | Hamming distances within and between clones |
| `*.starcode_.tsv` | Raw Starcode clustering results |
| `*.starcode_parsed.tsv` | Processed Starcode clustering output |
| `*.hX_corrected.umi_per_lineage.tsv` | Corrected UMIs per lineage barcode (after Starcode clustering) |
| `*.hX_corrected.lineages_per_cellbc.tsv` | Lineages per cell barcode after correction |
| `*.hX_corrected.clone_info.tsv` | Clone grouping after barcode correction |

🔹 **Legend for Naming Patterns:**
- `*` = `prefix` set by user (e.g., `Lin28_barcode`)
- `hX` = Maximum Hamming correction allowed (`--max_correction=X`)

## Troubleshooting
- **Missing Dependencies:** Ensure required Python packages are installed.
- **Starcode Not Found:** Verify that `--starcode_path` is correctly set.
- **Unexpected Output:** Check logs (`*.lin28_pipeline.log`) for errors.


