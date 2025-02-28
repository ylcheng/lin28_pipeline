import pandas as pd
import numpy as np
import argparse
import os
from itertools import combinations, product
# from scipy.spatial.distance import pdist, squareform
from datetime import datetime
import logging
import logging.config
import json


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
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def compute_distance_within(cellbc, lineages):
    """
    Compute Hamming distances within a single cellbc's lineage barcodes.

    Parameters:
    ----------
    cellbc : str
        The cell barcode.
    lineages : list of str
        List of lineage barcodes for this cellbc.

    Returns:
    -------
    list of dict
        Hamming distances within this cellbc.
    """
    distances = []
    for seq1, seq2 in combinations(lineages, 2):
        distances.append({
            "cellbc": cellbc,
            "lineage_bc1": seq1,
            "lineage_bc2": seq2,
            "hamming_distance": compute_hamming_distance(seq1, seq2),
            "comparison": "within_cellbc"
            })
    return distances

def compute_distance_between(cellbc, lineages, other_lineages):
    """
    Compute Hamming distances between this cellbc's lineages and other cellbcs' lineages

    Parameters:
    ----------
    cellbc : str
        The cell barcode.
    lineages : list of str
        List of lineage barcodes for this cellbc.
    other_lineages : list of str
        List of lineage barcodes from other cellbc.

    Returns:
    -------
    list of dict
        Hamming distances between this cellbc and other cellbcs.
    """
    distances = []
    for seq1, seq2 in product(lineages, other_lineages):
        distances.append({
            "cellbc": cellbc,
            "lineage_bc1": seq1,
            "lineage_bc2": seq2,
            "hamming_distance": compute_hamming_distance(seq1, seq2),
            "comparison": "between_cellbc",
            })
    return distances

def main(input_file, output_file,):
    """
    Main function to compute Hamming distances within and between cellbc lineage barcodes.
    """
    logger.info(f'Loading input file {input_file}...')
    df = pd.read_csv(input_file, sep = '\t')

    required_columns = ['cellbc','lineage_barcode' ]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        logger.error(f'Required columns not found in input file: {missing_columns}')
        raise ValueError(f'Required columns not found in input file: {missing_columns}')


    logger.info("Processing cellbcs for Hamming distance computation...")
    unique_cellbcs = df['cellbc'].unique()
    all_distances = []
    for index, cellbc in enumerate(unique_cellbcs, start=1):
        # Lineage barcodes for this cellbc
        lineages = df[df['cellbc'] == cellbc]['lineage_barcode'].tolist()
        # Lineage barcodes for all other cellbcs
        other_lineages = df[df['cellbc'] != cellbc]['lineage_barcode'].tolist()

        # Compute within-cellbc distances
        all_distances.extend(compute_distance_within(cellbc, lineages))
        # Compute between-cellbc distances
        all_distances.extend(compute_distance_between(cellbc, lineages, other_lineages))

        if index%1000 ==0:
            logger.info(f'Processed {index} cellbcs...')

    logger.info(f'Finished processing {len(unique_cellbcs)} cellbcs.')

    distances_df = pd.DataFrame(all_distances)
    # Save output in compressed format
    distances_df.to_csv(output_file, index = False, sep = '\t', compression = 'bz2')
    logger.info(f"Results saved  in compressed format to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Hamming distances of lineage barcodes within and between cellbc.")
    parser.add_argument("--input_file", required=True, help="Path to the input file (original umi_per_lineage or starcode corrected).")
    parser.add_argument("--output_file", required=True, help="Path to save the Hamming distance results. Save in Bzip2 format, so file extension should be bz2")
    parser.add_argument("--log_folder", default = None, help="Folder to save the log file.")
    parser.add_argument("--log_config", default = None, help="Path to the logging configuration file.")

    args = parser.parse_args()

    script_name = "compute_hamming_distances_within_and_between"
    logger = setup_logging( script_name, args.log_config, args.log_folder,)

    main(
        input_file=args.input_file,
        output_file=args.output_file,
    )






