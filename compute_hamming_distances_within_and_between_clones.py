import pandas as pd
import numpy as np
import argparse
import os
from itertools import combinations, product
from datetime import datetime
import logging
import logging.config
import json

def setup_logging(script_name, log_config=None, log_folder=None):
    """
    Configure logging using a JSON configuration file with a today_date log file name.
    """
    if log_config is None and log_folder is None:
        log_config = os.getenv('LOG_CONFIG')
        log_file = os.getenv("LOG_FILE")
    else:
        today_date = datetime.now().strftime("%Y-%m-%d")
        log_file = os.path.join(log_folder, f'{today_date}.{script_name}.log')
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

def compute_distance_within(cloneid, lineages):
    """
    Compute Hamming distances within a single clone's lineage barcodes.
    """
    distances = []
    for seq1, seq2 in combinations(lineages, 2):
        distances.append({
            "cloneid": cloneid,
            "lineage_bc1": seq1,
            "lineage_bc2": seq2,
            "hamming_distance": compute_hamming_distance(seq1, seq2),
            "comparison": "within_clone"
        })
    return distances

def compute_distance_between(cloneid, lineages, other_clone_lineages):
    """
    Compute Hamming distances between a clone's lineage barcodes and other clones' lineage barcodes.
    """
    distances = []
    for seq1, seq2 in product(lineages, other_clone_lineages):
        distances.append({
            "cloneid": cloneid,
            "lineage_bc1": seq1,
            "lineage_bc2": seq2,
            "hamming_distance": compute_hamming_distance(seq1, seq2),
            "comparison": "between_clones"
        })
    return distances

def main(input_file, output_file):
    """
    Main function to compute Hamming distances within and between cloneid lineage barcodes.
    """
    logger.info(f'Loading input file {input_file}...')
    df = pd.read_csv(input_file, sep='\t')

    required_columns = ['cloneid', 'sorted_lineage_barcodes']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        logger.error(f"Required columns not found in input file: {missing_columns}")
        raise ValueError(f"Required columns not found in input file: {missing_columns}")

    # Preprocess sorted_lineage_barcodes to create individual rows for each barcode
    logger.info("Expanding sorted_lineage_barcodes into individual rows for processing...")
    df['lineage_barcodes_list'] = df['sorted_lineage_barcodes'].apply(lambda x: x.split(','))
    df = df.explode('lineage_barcodes_list').rename(columns={'lineage_barcodes_list': 'lineage_barcode'})

    logger.info("Processing clones for Hamming distance computation...")
    unique_cloneids = df['cloneid'].unique()
    all_distances = []

    for index, cloneid in enumerate(unique_cloneids, start=1):
        # Lineage barcodes for this cloneid
        # lineages = df[df['cloneid'] == cloneid]['sorted_lineage_barcodes'].iloc[0].split(',')
        lineages = df[df['cloneid'] == cloneid]['lineage_barcode'].tolist()
        # Lineage barcodes for all other cloneids
        # other_lineages = df[df['cloneid'] != cloneid]['sorted_lineage_barcodes'].str.split(',').explode().tolist()
        other_lineages = df[df['cloneid'] != cloneid]['lineage_barcode'].tolist()

        # Compute within-clone distances
        all_distances.extend(compute_distance_within(cloneid, lineages))
        # Compute between-clone distances
        all_distances.extend(compute_distance_between(cloneid, lineages, other_lineages))

        if index % 1000 == 0:
            logger.info(f'Processed {index} clones...')

    logger.info(f'Finished processing {len(unique_cloneids)} clones.')

    distances_df = pd.DataFrame(all_distances)
    # Save output in compressed format
    distances_df.to_csv(output_file, index=False, sep='\t', compression='bz2')
    logger.info(f"Results saved in compressed format to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Hamming distances of lineage barcodes within and between clones.")
    parser.add_argument("--input_file", required=True, help="Path to the input file (output from create_clone_info).")
    parser.add_argument("--output_file", required=True, help="Path to save the Hamming distance results. Save in Bzip2 format, so file extension should be bz2.")
    parser.add_argument("--log_folder", default=None, help="Folder to save the log file.")
    parser.add_argument("--log_config", default=None, help="Path to the logging configuration file.")

    args = parser.parse_args()

    script_name = "compute_hamming_distances_within_and_between_clones"
    logger = setup_logging(script_name, args.log_config, args.log_folder)

    main(
        input_file=args.input_file,
        output_file=args.output_file
    )

