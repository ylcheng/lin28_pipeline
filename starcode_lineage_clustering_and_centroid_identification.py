import pandas as pd
import os
import subprocess
import json
import logging
import logging.config
from datetime import datetime
import argparse
from itertools import product
import numpy as np



def setup_logging( script_name, log_config=None, log_folder=None,):
    """
    Configure logging using a JSON configuration file with a today_date log file name.
    """
    if log_config is None and log_folder is None:
        log_config = os.getenv('LOG_CONFIG')
        log_file = os.getenv("LOG_FILE")
    else:
        today_date = datetime.now().strftime("%Y-%m-%d")
        # today_date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        log_file = os.path.join(log_folder, f'{today_date}.{script_name}.log')
    print(log_config)
    print(log_folder)
    with open(log_config, 'r') as config_file:
        config = json.load(config_file)
    config['handlers']['file']['filename'] = log_file
    logging.config.dictConfig(config)
    return logging.getLogger()



def compute_hamming_distance(seq1, seq2):
    """
    Compute Hamming distance between two sequences.
    """
    if len(seq1) != len(seq2):
        return np.nan
    else:
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def run_starcode(input_file, output_file, starcode_path, distance, cluster_method , cluster_ratio , index):
    """
    Run Starcode on lineage barcodes in each cell barcodes

    --distance 
     Defines the maximum Levenshtein distance for clustering.  When not set it is automatically computed as: min(8, 2 + [median seq length]/30)

     --cluster-ratio
    (Message passing only) Specifies the minimum sequence count ratio to cluster two matching sequences, i.e. two matching sequences A and B will be clustered together only if count(A) > ratio * count(B).
     Sparse datasets may need to set -r to small values (minimum is 1.0) to trigger clustering.
     Default is 5.0.

     --sphere
      Use sphere clustering algorithm instead of message passing (MP). Spheres is more greedy than MP: sorted by size, centroids absorb all their matches.

    --print-clusters
    Adds a third column to the starcode output, containing the sequences that compose each cluster.
    By default, the output contains only the centroid and the counts.


    """
    # baseline command
    command = [
            starcode_path,
            "--input", input_file,
            "--output", output_file,
            "--dist", str(distance),
            "--print-clusters",
            ]
    if cluster_method == 'sphere':
        command = command + ["--sphere"]
    elif cluster_method == 'message_passing':
        command = command + ["--cluster-ratio", cluster_ratio]

    # Run starcode on terminal
    # Log first and every 1000 cellls
    if index == 1 or index%1000 ==0:
        logger.info(f'Running Starcode with command: {" ".join(command)}')
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    if index == 1 or index%1000 ==0:
        logger.info(f"STDOUT: {result.stdout}")
        logger.error(f"STDERR: {result.stderr}")
        logger.info(f'Starcode output saved to {output_file}') 


def process_starcode_output(cellbc, starcode_output_file, compute_distance ):
    """
    Parse Starcode output and generate DataFrame rows for each member of each cluster.
    """
    rows = []
    with open(starcode_output_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            centroid = parts[0] # canonical sequence
            barcodes = parts[2].split(',') # cluster sequences
            for barcode in barcodes:
                hamming_distance = compute_hamming_distance(barcode, centroid) if compute_distance else None
                rows.append({
                    "cellbc": cellbc,
                    "lineage_barcode": barcode,
                    "centroid": centroid,
                    "hamming_distance": hamming_distance
                    })
    # logger.info(f'Processed {len(rows)} rows for cellbc {cellbc}') # testing codes, comment out later
    return rows

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
    if cluster_method == "messagePassing":
        # Add cluster ratio for message passing method
        suffix += f"_r{cluster_ratio}"

    # Replace the file extension and append the suffix
    # base, ext = os.path.splitext(output_file)
    # return f"{base}{suffix}{ext}"
    return output_file.replace(".tsv", f'{suffix}.tsv')



def main(umi_per_lineage_file, output_file, cluster_method, distance, cluster_ratio, compute_distance,
         folder_path, starcode_path, cleanup_temp, combined_starcode_output):
    """
    Main function for clustering lineage barcodes per cell barcode using Starcode.

    This function performs the following steps:
    1. Reads an input file containing cell barcodes (`cellbc`), lineage barcodes, and their UMI counts.
    2. Expands lineage barcodes by their UMI counts and prepares temporary input files for each `cellbc`.
    3. Runs Starcode clustering on the lineage barcodes of each `cellbc`.
    4. Parses the Starcode output to create a structured DataFrame containing:
       - `cellbc`: The cell barcode.
       - `lineage_barcode`: A barcode in the cluster.
       - `centroid`: The canonical sequence (cluster centroid).
       - `hamming_distance`: Optional, distance between the barcode and its centroid.
    5. Optionally collects all Starcode outputs into a combined file, adding a `cellbc` column.
    6. Cleans up temporary files.

    Parameters:
    ----------
    umi_per_lineage_file : str
        Path to the input file containing columns `cellbc`, `lineage_barcode`, and `umi_counts`.
    output_file : str
        Path to save the final output with parsed Starcode results.
    cluster_method : str
        Clustering method for Starcode. Choices: `sphere` or `message_passing`.
    distance : int
        Maximum Levenshtein distance allowed for clustering.
    cluster_ratio : float
        Minimum sequence count ratio for merging clusters in `message_passing` clustering. Ignored if `sphere` is selected.
    compute_distance : bool
        If True, computes and includes the Hamming distance between lineage barcodes and their centroids.
    folder_path : str
        Directory to save temporary Starcode output per cellbc files and the combined Starcode output file.
    starcode_path : str
        Path to the Starcode executable.
    cleanup_temp : bool
        If True, deletes temporary files and directories created during processing.
    combined_starcode_output : str or None
        Name of the combined Starcode output file. If None, this feature is disabled.
    """

    

    logger.info(f'Loading umi_per_lineage_file {umi_per_lineage_file}...')
    umi_df = pd.read_csv(umi_per_lineage_file, sep = '\t')

    final_output = []
    combined_starcode_rows = ['cellbc\tcentroid\tsize\tsequences\n']
    combined_starcode_rows = []
    temp_dir = f'{folder_path}/temp_starcode_files'
    os.makedirs(temp_dir, exist_ok = True)
    logger.info(f'Temporary folder for starcode output files: {temp_dir}...')

    unique_cellbcs = umi_df['cellbc'].unique() # uncomment 
    # unique_cellbcs = umi_df['cellbc'].unique()[:3] # testing codes, comment out later
    for index, cellbc in enumerate(unique_cellbcs, start = 1):
        cellbc_df = umi_df[umi_df['cellbc'] == cellbc]
        lineage_barcodes = []
        # lineage barcode represented by its umi_counts
        for _, row in cellbc_df.iterrows():
            lineage_barcodes.extend([row['lineage_barcode']] * row['umi_counts'])

        temp_input_file = os.path.join(temp_dir, f'{cellbc}_input.txt')
        temp_output_file = os.path.join(temp_dir, f'{cellbc}_output.txt')

        # write temporary Starcode input file for lineage barcodes of the given cellbc
        with open(temp_input_file, 'w') as file:
            file.write("\n".join(lineage_barcodes))

        # Run Starcode on the lineage barocdes of the given cellbc
        run_starcode(temp_input_file, temp_output_file, starcode_path, distance, cluster_method , cluster_ratio , index)
        # Process Starcode output for the final parsed structure
        final_output.extend(process_starcode_output(cellbc, temp_output_file, compute_distance))

        # Optionally collect rows for the combined Starcode output
        if combined_starcode_output:
            with open(temp_output_file, 'r') as file:
                for line in file:
                    parts = line.strip().split("\t")
                    combined_starcode_rows.append({
                        "cellbc": cellbc,
                        "centroid" : parts[0],
                        "size" : parts[1],
                        "sequences" : parts[2],
                        })

        # Log progress statement every 1000 cells
        if index % 1000 ==0:
            logger.info(f'Processed {index} cells...')
    logger.info(f'Finished starcode and parsing; Processed {len(unique_cellbcs)} cells')

   # Add clustering parameters to suffix of output_file
    # output_file = construct_output_filename(output_file, cluster_method, distance, cluster_ratio)

    # Convert final output to Dataframe and save
    final_df = pd.DataFrame(final_output)
    final_df.to_csv(output_file, sep = '\t', index = False)
    logger.info(f'Final parsed Starcode output saved {output_file}')

    # Save combined Starcode output if requested
    if combined_starcode_output:
        # combined_starcode_output = construct_output_filename(combined_starcode_output, cluster_method, distance, cluster_ratio)
        combined_df = pd.DataFrame(combined_starcode_rows)
        combined_df.to_csv(combined_starcode_output, sep = '\t', index = False)
        logger.info(f'Combined Starcode output saved to {combined_starcode_output}')

    # Cleanup termporary Starcode files
    if cleanup_temp:
        for temp_file in os.listdir(temp_dir):
            os.remove(os.path.join(temp_dir, temp_file))
        os.rmdir(temp_dir)
        logger.info("Starcode temporary files cleaned up.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster lineage barcodes using Starcode for each cell barcode and identify canonical sequences (centroids).")
    # Input and Output Files
    parser.add_argument("--umi_per_lineage_file", required=True, help="Path to the input file containing columns 'cellbc', 'lineage_barcode', and 'umi_counts'.")
    parser.add_argument("--output_file", required=True, help="Path to save the final output with clustered lineage barcodes and their centroids.")
    parser.add_argument("--combined_starcode_output", default=None, help="Path to save the combined Starcode output with cellbc column. If not set, combined output is not generated.")
    # Logging Configuration
    # Clustering Parameters
    parser.add_argument("--cluster_method", choices=["sphere", "message_passing"], default="sphere", help="Clustering method to use with Starcode. Default is 'sphere'.")
    parser.add_argument("--distance", type=int, default = 1,  help="Maximum Levenshtein distance allowed for clustering sequences in Starcode.")
    parser.add_argument("--cluster_ratio", type=float, default=5.0, help="Minimum sequence count ratio for merging clusters in 'message_passing' clustering. Ignored if 'sphere' clustering is selected. Default is 5.0.")
    # Additional Processing Options
    parser.add_argument("--compute_distance", action="store_true", help="If set, computes and includes the Hamming distance between lineage barcodes and their centroids in the final output.")
    parser.add_argument("--folder_path", required=True, help="Directory to save temporary input and output files generated during processing. Temporary files will be cleaned up if the '--cleanup_temp' flag is used.")
    parser.add_argument("--starcode_path", required=True, help="Path to the Starcode executable. Ensure Starcode is installed and accessible.")
    parser.add_argument("--cleanup_temp", action="store_true", help="If set, deletes all temporary files and directories created during processing.")
    parser.add_argument("--log_folder", default = None, help="Folder to save the log file.")
    parser.add_argument("--log_config", default = None, help="Path to the logging configuration file.")

    args = parser.parse_args()

    script_name = "starcode_lineage_clustering_and_centroid_identification" # no file extension .py
    logger = setup_logging( script_name, args.log_config, args.log_folder,)

    # Run the main function with the parsed arguments
    main(
        umi_per_lineage_file=args.umi_per_lineage_file,
        output_file=args.output_file,
        cluster_method=args.cluster_method,
        distance=args.distance,
        cluster_ratio=args.cluster_ratio,
        compute_distance=args.compute_distance,
        folder_path=args.folder_path,
        starcode_path=args.starcode_path,
        cleanup_temp=args.cleanup_temp,
        combined_starcode_output = args.combined_starcode_output,
    )

