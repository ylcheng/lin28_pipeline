import os
import subprocess
import logging
import logging.config
import json
from datetime import datetime
import argparse



def run_subprocess(command, script):
    """
    Run a subprocess and log its progress.
    """
    try:
        start_time = datetime.now()
        logger.info(f"EXECUTING: {' '.join(command)}")
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        logger.info(f"STDOUT: {result.stdout}")
        logger.error(f"STDERR: {result.stderr}")
        end_time = datetime.now()
        logger.info(f'COMPLETED {script} in {end_time - start_time}')
    except subprocess.CalledProcessError as error:
        logger.error(f'FAILED {script}: {error}')
        logger.error(f"STDERR: {error.stderr}")
        raise


def lin28_pipeline(args):
    """
    Run the Lin28 lineage barcode analysis pipeline to process sequencing data,
    compute lineage barcodes, and analyze clonal relationships.

    Steps:
    -------
    1. Subset Reads with Flanking Sequences:
       Identify and extract reads containing specific start and stop flanking sequences in Read 2.
       Save the corresponding paired reads (Read 1 and Read 2) to new FASTQ files.

    2. Extract Lineage Barcode Data:
       Parse the extracted FASTQ files to identify cell barcodes, unique molecular identifiers (UMIs),
       and lineage barcodes. Additionally, compute mismatches for each sequence relative to the expected sequences.

    3. Compute Hamming Distance Between Cell Barcodes (if max_hamming_distance > 0):
       Calculate Hamming distances between cell barcodes in the RNA dataset and the lineage dataset.
       Enables either all-pairwise or segment-based matches, depending on user settings.

    4. Filter Lineage Barcode Reads and Correct Cell Barcodes:
       Apply thresholds for mismatches in cell barcodes and lineage barcodes.
       Filter out reads and correct cell barcodes based on RNA dataset references and max_hamming_distance.

    5. Collapse Reads Per UMI:
       Combine reads with the same cell barcode, UMI, and lineage barcode.
       Count the total reads per UMI.

    6. Collapse UMIs Per Lineage Barcode:
       For each cell barcode, count the number of UMIs associated with each lineage barcode.
       Count the number of UMIs per lineage barcode.

    7. Generate Lineages Per Cell Barcode (if generate_clones_associated_files_before_correction is True or max_correction = 0):
       Create the summary file containing the number of unique lineage barcodes and their sequences
       for each cell barcode.

    8. Assign Clone IDs to Unique Lineage Barcode Sets (if generate_clones_associated_files_before_correction is True or max_correction = 0):
       Group cells by their unique sets of lineage barcodes and assign a clone ID.
       Record the number of cells and lineage barcodes for each clone.

    9. Compute Hamming Distances Within and Between Clones (if compute_hamming_distances_within_and_between_clones is True):
       Analyze Hamming distances between lineage barcodes within individual clones
       and across different clones to assess sequence similarity.

    10. Cluster Lineage Barcodes Using Starcode (if max_correction > 0):
        Perform clustering on lineage barcodes using Starcode. Identify centroids and their associated
        lineage barcodes that are at most levenshtein_distance away.

    11. Correct and Collapse UMIs Per Lineage Barcode (if max_correction > 0):
        Use the Starcode clustering results to correct lineage barcodes and recalculate
        UMIs for each lineage barcode.

    12. Generate Lineages Per Cell Barcode After Correction (if max_correction > 0):
        Create the summary file of lineage barcodes per cell barcode after applying Starcode corrections.

    13. Assign Clone IDs to Corrected Lineage Barcode Sets (if max_correction > 0):
        Group cells by unique sets of corrected lineage barcodes and assign a clone ID.
        Record the number of cells and lineage barcodes for each clone.


    Parameters:
    -----------
    args : argparse.Namespace
        Parsed command-line arguments containing input/output paths, thresholds, and configurations.

    Returns:
    --------
    None
    """
    logger.info(f"Executing pipeline from step {args.start_step} to step {args.end_step}")

    # STEP 1: Subset reads containing start and stop flanking sequences in Read 2.
    # Save paired reads (Read 1 and Read 2) with detected flanking sequences to new FASTQ files.
    if args.start_step <= 1 <= args.end_step:
        script = "subset_reads_with_flanking_sequences.py"
        command = [
                "python", f"{args.code_path}/{script}",
                "--folder_path", args.folder_path,
                "--output_fastq_prefix", args.prefix,
                "--start_sequence", args.start_flanking,
                "--stop_sequence", args.stop_flanking,
                "--fastq_pattern", args.fastq_pattern,
                "--R1_pattern", args.R1_pattern,
                "--R2_pattern", args.R2_pattern,
                ]

        run_subprocess(command = command, script = script)


    # STEP 2: Extract lineage barcodes, cell barcodes, and UMIs from the subset FASTQ files.
    # Compute mismatches for each sequence relative to expected flanking and barcode sequences.
    # output: read_id,	cellbc,	umi,	lineage_barcode,	number_of_mismatch_lineage_barcode,	start_flanking,	number_of_mismatch_start_flanking,	stop_flanking,	number_of_mismatch_stop_flanking
    path_prefix = os.path.join(args.folder_path, args.prefix)
    lineage_reads_file = f"{path_prefix}.lineage_barcode_reads.tsv"
    troubleshoot_file = f"{path_prefix}.ambiguous_lineage_barcode_reads.tsv"
    r1_file = f"{path_prefix}_R1_flanking_detected.fastq.gz"
    r2_file = f"{path_prefix}_R2_flanking_detected.fastq.gz"

    if args.start_step <= 2 <= args.end_step:
        script = "extract_lineage_barcode_data.py"
        command = [
            "python", f"{args.code_path}/{script}",
            "--r1_file", r1_file,
            "--r2_file", r2_file,
            "--output_file", lineage_reads_file,
            "--troubleshoot_file", troubleshoot_file,
            "--start_flanking", args.start_flanking,
            "--stop_flanking", args.stop_flanking
        ]
        run_subprocess(command = command, script = script)


    # STEP 3: Compute Hamming distances between cell barcodes of RNA dataset and lineage dataset
    # Perform either all-pairwise or segment-based pre-filter before computation
    # output: Reference,	Query,	HammingDistance,
    hamming_file = f'{args.folder_path}/{args.prefix}.cell_barcodes_hamming_distances.tsv.bz2'

    if args.start_step <= 3 <= args.end_step:
        script = "compute_hamming_distances_btw_references_and_queries.py"

        if args.max_hamming_distance > 0 and args.segment_match:
            # only compute those combinations with potential Hamming distance being at most max_hamming_distance
            # only kept those combinations with Hamming distance less than or equal to max_hamming_distance
            logger.info("Segment match filtering enabled.")
            command = [
                    "python", f'{args.code_path}/{script}',
                    "--reference_file",  args.rna_file,
                    "--query_file",  lineage_reads_file,
                    "--output_file",  hamming_file,
                    "--max_hamming_distance", str(args.max_hamming_distance),
                    "--segment_match",
                    ]
            run_subprocess(command = command, script = script)
        elif args.all_pairwise:
            # computing and kept all pairwise combinations
            logger.info("All pairwise computation enabled.")
            command = [
                    "python", f'{args.code_path}/{script}',
                    "--reference_file",  args.rna_file,
                    "--query_file",  lineage_reads_file,
                    "--output_file",  hamming_file,
                    ]

            run_subprocess(command = command, script = script)


    # STEP 4: Apply mismatch thresholds to filter lineage barcodes and correct cell barcodes using RNA cell barcodes as reference
    if args.max_hamming_distance > 0 and args.segment_match:
        hamming_file = hamming_file.replace(".tsv", f".max_hamming_{args.max_hamming_distance}.tsv")
    elif args.all_pairwise:
        hamming_file = hamming_file.replace(".tsv", ".all_pairwise.tsv")
    else:
        hamming_file = None

    filtered_reads_file = f'{args.folder_path}/{args.prefix}.filtered_lineage_barcode_reads.tsv'

    if args.start_step <= 4 <= args.end_step:
        script = "filter_lineage_barcode_reads.py"
        if hamming_file is None:
            command = [
                    "python", f'{args.code_path}/{script}',
                    "--valid_file", lineage_reads_file,
                    "--rna_file", args.rna_file,
                    "--max_hamming_distance", str(args.max_hamming_distance),
                    "--max_lineage_mismatch", str(args.max_lineage_mismatch),
                    "--max_start_flanking_mismatch", str(args.max_start_flanking_mismatch),
                    "--max_stop_flanking_mismatch", str(args.max_stop_flanking_mismatch),
                    "--output_file", filtered_reads_file,
                    ]
            run_subprocess(command = command, script = script)
        else:
            command = [
                    "python", f'{args.code_path}/{script}',
                    "--valid_file", lineage_reads_file,
                    "--hamming_file", hamming_file,
                    "--rna_file", args.rna_file,
                    "--max_hamming_distance", str(args.max_hamming_distance),
                    "--max_lineage_mismatch", str(args.max_lineage_mismatch),
                    "--max_start_flanking_mismatch", str(args.max_start_flanking_mismatch),
                    "--max_stop_flanking_mismatch", str(args.max_stop_flanking_mismatch),
                    "--output_file", filtered_reads_file,
                    ]
            run_subprocess(command = command, script = script)



    # STEP 5: Combine reads with the same cell barcode, UMI, and lineage barcode.
    # Count the total number of reads for each UMI.
    # output: cellbc,	umi,	lineage_barcode,	reads_count
    output_suffix = (f".cellbc{args.max_hamming_distance}"
                 f"start{args.max_start_flanking_mismatch}"
                 f"lineage{args.max_lineage_mismatch}"
                 f"stop{args.max_stop_flanking_mismatch}")
    filtered_reads_file = filtered_reads_file.replace(".tsv", f"{output_suffix}.tsv")
    reads_per_umi_file = filtered_reads_file.replace(".tsv", f".reads_per_umi.tsv")

    if args.start_step <= 5 <= args.end_step:
        script = "collapse_reads_per_umi.py"
        command = [
            "python", f"{args.code_path}/{script}",
            "--input_file", filtered_reads_file,
            "--output_file", reads_per_umi_file
        ]
        run_subprocess(command = command, script = script)


        # STEP 6: Collapse unique UMIs for each cell barcode and lineage barcode combination.
        # Calculate the total number of UMIs per lineage barcode.
        # output: cellbc,	lineage_barcode,	umi_counts
        umi_per_lineage_file = filtered_reads_file.replace(".tsv", ".umi_per_lineage.tsv")
        script = "collapse_umi_per_lineage.py"
        command = [
            "python", f"{args.code_path}/{script}",
            "--input_file", reads_per_umi_file,
            "--output_file", umi_per_lineage_file
        ]
        run_subprocess(command = command, script = script)


    if args.generate_clones_associated_files_before_correction is True or args.max_correction == 0:
        # STEP 7: Summarize the number and sequences of unique lineage barcodes detected for each cell barcode.
        # output: cellbc,	lineages_count,	lineage_barcodes,	sorted_lineage_barcodes,	cloneid(append after STEP 8)
        lineages_per_cellbc_file = filtered_reads_file.replace(".tsv", ".lineages_per_cellbc.tsv")

        if args.start_step <= 7 <= args.end_step:
            script = "generate_lineages_per_cellbc.py"
            command = [
                "python", f"{args.code_path}/{script}",
                "--input_file", umi_per_lineage_file,
                "--output_file", lineages_per_cellbc_file
            ]
            run_subprocess(command = command, script = script)

        # STEP 8: Group clones by unique sets of lineage barcodes and assign clone IDs.
        # Record the number of cells and lineage barcodes for each clone.
        # output: sorted_lineage_barcodes,	lineages_count,	number_of_cells,	cloneid,
        clone_info_file = filtered_reads_file.replace('.tsv', f'.clone_info.tsv')

        if args.start_step <= 8 <= args.end_step:
            script = "create_clone_info.py"
            command = [
                    "python", f'{args.code_path}/{script}',
                    "--input_file", lineages_per_cellbc_file,
                    "--output_file", clone_info_file,
                    ]
            run_subprocess(command = command, script = script)


    if args.compute_hamming_distances_within_and_between_clones is True:
        # STEP 9: Analyze Hamming distances within individual clones and across different clones.
        # Assess sequence similarity and overlap among lineage barcodes between clones
        # output: cloneid,	lineage_bc1,	lineage_bc2,	hamming_distance,	comparison
        clone_info_file = filtered_reads_file.replace('.tsv', f'.clone_info.tsv')
        lineage_hamming_file = filtered_reads_file.replace(".tsv", f'.lineage_hammings_within_and_btw_clones.tsv.bz2')
        if args.start_step <= 9 <= args.end_step:
            if not os.path.isfile(clone_info_file):
                logger.error(f"Not computing Hamming distance of lineage barcodes within and between clones, because input file does not exist: {clone_info_file}")
                raise FileNotFoundError(f"Not computing Hamming distance of lineage barcodes within and between clones, because input file does not exist: {clone_info_file}")
            script = "compute_hamming_distances_within_and_between_clones.py"
            command = [
                    "python", f'{args.code_path}/{script}',
                    "--input_file", clone_info_file,
                    "--output_file", lineage_hamming_file,
                    ]
            run_subprocess(command = command, script = script)


    if args.max_correction > 0:
        # STEP 10: Perform clustering of lineage barcodes using Starcode.
        # Identify centroids and their associated lineage barcodes that are at most levenshtein_distance away.
        # output: cellbc,	lineage_barcode,	centroid,	hamming_distance
        def construct_output_filename(output_file, cluster_method, distance, cluster_ratio=None):
            """
            Construct the output filename with a concise suffix based on clustering parameters.

            Parameters:
            ----------
            output_file : str
                Base output file name.
            cluster_method : str
                Clustering method (e.g., "sphere" or "message_passing").
            distance : int
                Maximum Levenshtein distance.
            cluster_ratio : float
                Cluster ratio (used only for "message_passing").

            Returns:
            -------
            str
                Modified output file name with suffix.
            """
            # Base suffix for cluster method and distance
            suffix = f"{cluster_method}_d{distance}"
            if cluster_method == "message_passing":
                # Add cluster ratio for message passing method
                suffix += f"_r{cluster_ratio}"

            # Replace the file extension and append the suffix
            return output_file.replace(".tsv", f'{suffix}.tsv')

        starcode_file = filtered_reads_file.replace(".tsv", ".starcode_.tsv")
        starcode_file = construct_output_filename(starcode_file, args.cluster_method, args.levenshtein_distance)
        parsed_starcode_file = starcode_file.replace(".tsv", "_parsed.tsv")

        if args.start_step <= 10 <= args.end_step:
            script = "starcode_lineage_clustering_and_centroid_identification.py"
            command = [
                    "python", f'{args.code_path}/{script}',
                    "--umi_per_lineage_file", umi_per_lineage_file,
                    "--output_file", parsed_starcode_file,
                    "--combined_starcode_output", starcode_file,
                    "--distance", str(args.levenshtein_distance),
                    "--compute_distance",
                    "--folder_path", args.folder_path,
                    "--starcode_path", args.starcode_path,
                    "--cleanup_temp",
                    ]
            run_subprocess(command = command, script = script)


        # STEP 11: Use Starcode clustering results to correct lineage barcodes.
        # Collapse UMI counts for the corrected lineage barcodes.
        # output: cellbc,	lineage_barcode,	umi_counts
        starcode_corrected_umi_per_lineage_file = starcode_file.replace('.tsv', f'.h{args.max_correction}_corrected.umi_per_lineage.tsv')
        if args.start_step <= 11 <= args.end_step:
            script = "correct_lineages_with_starcode_and_collapse_umi.py"
            command = [
                    "python", f'{args.code_path}/{script}',
                    "--centroid_file", parsed_starcode_file,
                    "--umi_file", umi_per_lineage_file,
                    "--output_file", starcode_corrected_umi_per_lineage_file,
                    "--max_hamming_distance", str(args.max_correction),
                    ]
            run_subprocess(command = command, script = script)


        # STEP 12: Summarize the number and sequences of corrected lineage barcodes for each cell barcode.
        # output: cellbc,	lineages_count,	lineage_barcodes,	sorted_lineage_barcodes,	cloneid(append after STEP 13)
        starcode_corrected_lineages_per_cellbc_file = starcode_file.replace('.tsv', f'.h{args.max_correction}_corrected.lineages_per_cellbc.tsv')

        if args.start_step <= 12 <= args.end_step:
            script = "generate_lineages_per_cellbc.py"
            command = [
                    "python", f'{args.code_path}/{script}',
                    "--input_file", starcode_corrected_umi_per_lineage_file,
                    "--output_file", starcode_corrected_lineages_per_cellbc_file,
                    ]
            run_subprocess(command = command, script = script)


        # STEP 13: Group cells by unique sets of corrected lineage barcodes and assign a clone ID
        # Record the number of cells and lineage barcodes for each corrected clone.
        # output: sorted_lineage_barcodes,	lineages_count,	number_of_cells,	cloneid
        starcode_corrected_clone_info_file = starcode_file.replace('.tsv', f'.h{args.max_correction}_corrected.clone_info.tsv')
        if args.start_step <= 13 <= args.end_step:
            script = "create_clone_info.py"
            command = [
                    "python", f'{args.code_path}/{script}',
                    "--input_file", starcode_corrected_lineages_per_cellbc_file,
                    "--output_file", starcode_corrected_clone_info_file,
                    ]
            run_subprocess(command = command, script = script)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the Lin28 lineage barcode analysis pipeline to process sequencing data, compute lineage barcodes, and analyze clonal relationships.")

    # General paths and file patterns
    parser.add_argument("--code_path", default=".", help="Path to the directory containing the pipeline scripts (default: current directory).")
    parser.add_argument("--folder_path", required=True, help="Path to the input/output directory for processing.")
    parser.add_argument("--prefix", required=True, help="Base prefix for naming output files.")
    parser.add_argument("--fastq_pattern", default="*fastq.gz", help="Filename pattern to search for gzipped FASTQ files in the folder (default: '*fastq.gz').")
    parser.add_argument("--R1_pattern", default="_R1_", help="Pattern to identify Read 1 (R1) files in FASTQ files (default: '_R1_').")
    parser.add_argument("--R2_pattern", default="_R2_", help="Pattern to identify Read 2 (R2) files in FASTQ files (default: '_R2_').")

    # Flanking sequences
    parser.add_argument("--start_flanking", default="GAAACACCG", help="Sequence to identify the start flanking region (default: 'GAAACACCG').")
    parser.add_argument("--stop_flanking", default="GTTTTAGAG", help="Sequence to identify the stop flanking region (default: 'GTTTTAGAG').")

    # Mismatch thresholds
    parser.add_argument("--max_hamming_distance", type=int, default=0, help="Maximum allowable Hamming distance for matching cell barcodes (default: 0, exact match). Based on 10x Genomics whitelist, can set to 1")
    parser.add_argument("--max_lineage_mismatch", type=int, default=0, help="Maximum allowable mismatches for lineage barcodes (default: 0, exact match).")
    parser.add_argument("--max_start_flanking_mismatch", type=int, default=1, help="Maximum allowable mismatches in the start flanking sequence (default: 1).")
    parser.add_argument("--max_stop_flanking_mismatch", type=int, default=1, help="Maximum allowable mismatches in the stop flanking sequence (default: 1).")

    # RNA dataset reference
    parser.add_argument("--rna_file", required=True, help="Path to the RNA dataset file containing reference cell barcodes.")
    # parser.add_argument("--reference_column", type=int, default=0, help="Index of the column (0-based) in the RNA dataset containing reference cell barcodes (default: 0).")
    # parser.add_argument("--reference_separator", default="\t", help="Column separator used in the RNA dataset file (default: tab).")
    # parser.add_argument("--reference_header", default=None, help="Specify if the RNA dataset file contains a header. Use 'None' for no header (default: None).")

    # Cell barcode Hamming distance computation
    parser.add_argument("--segment_match", action="store_true", help="Enable segment match filtering to identify potential cell barcode matches based on Hamming distance thresholds.")
    parser.add_argument("--all_pairwise", action="store_true", help="Enable all-pairwise comparison for cell barcode matching.")

    # Clone-related processing before lineage barcode correction for comparision
    parser.add_argument("--generate_clones_associated_files_before_correction", action="store_true", help="Generate lineage_per_cellbc and clone information files before lineage barcode correction.")
    parser.add_argument("--compute_hamming_distances_within_and_between_clones", action="store_true", default=False, help="Enable computation of Hamming distances within and between clones (default: False).")

    # Starcode clustering for lineage barcode correction
    parser.add_argument("--starcode_path", required=True, help="Path to the Starcode executable. Ensure Starcode is installed and accessible.")
    parser.add_argument("--cluster_method", choices=["sphere", "message_passing"], default="sphere", help="Clustering method to use in Starcode: 'sphere' for distance-based clustering or 'message_passing' for iterative merging (default: 'sphere').")
    parser.add_argument("--levenshtein_distance", type=int, default=2, help="Maximum Levenshtein distance for clustering sequences in Starcode (default: 2).")
    parser.add_argument("--cluster_ratio", type=float, default=5.0, help="Minimum sequence count ratio for merging clusters in 'message_passing' mode (ignored for 'sphere' clustering; default: 5.0).")
    parser.add_argument("--max_correction", type=int, default=2, help="Maximum number of nucleotide corrections for lineage barcodes; (estimated minimum Hamming distance - 1) / 2 (default: 2,  Based on minimum Hamming distance of lineages barcodes between clones in one observed data is about 5")

    # Control STEP to execute
    parser.add_argument("--start_step", type=int, default=1, help="Step number to start the pipeline from (default: 1).")
    parser.add_argument("--end_step", type=int, default=13, help="Step number to stop the pipeline at (default: 13).")

    args = parser.parse_args()

    # Validate mutually exclusive arguments and logical conditions for STEP 3: Compute cell barcodes Hamming distance 
    if args.max_hamming_distance == 0:
        if args.segment_match or args.all_pairwise:
            parser.error("With max_hamming_distance = 0, neither --segment_match nor --all_pairwise should be enabled.")
    elif args.max_hamming_distance > 0:
        if args.segment_match and args.all_pairwise:
            parser.error("Do not enable both --segment_match and --all_pairwise simultaneously when max_hamming_distance > 0.")
        if not (args.segment_match or args.all_pairwise):
            parser.error("With max_hamming_distance > 0, you must enable either --segment_match or --all_pairwise.")


    # STEP 0: Setup Logging
    # Get today's date in the desired format
    # today_date = datetime.now().strftime("%Y-%m-%d")
    today_date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file = f"{args.folder_path}/{today_date}.lin28_pipeline.log"
    log_config = f"{args.code_path}/logging_config.json"
    # Set environment variables for logging
    os.environ["LOG_CONFIG"] = log_config
    os.environ["LOG_FILE"] = log_file
    # Load the logging configuration in the parent
    with open(os.environ["LOG_CONFIG"], "r") as config_file:
        config = json.load(config_file)
    config["handlers"]["file"]["filename"] = os.environ["LOG_FILE"]
    config["handlers"]["file"]["mode"] = "w"
    logging.config.dictConfig(config)

    logger = logging.getLogger()

    pipeline_start_time = datetime.now()
    logger.info('START lin28_pipeline')
    logger.info(f'LOG_FILE: {log_file}')
    logger.info(f'LOG_CONFIG: {log_config}')
    # Log arguments
    logger.info('Pipeline Arguments:')
    for arg, value in vars(args).items():
        logger.info(f"{arg}: {value}")


    lin28_pipeline(args)

    pipeline_end_time = datetime.now()
    logger.info('END lin28_pipeline')
    logger.info(f'Total Run Time: {pipeline_end_time - pipeline_start_time}')

