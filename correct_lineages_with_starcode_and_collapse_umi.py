import pandas as pd
import argparse
import logging
import logging.config
import os
import json
from datetime import datetime


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


def correct_and_collapse_umi_per_lineage(centroid_file, umi_file, output_file, max_hamming_distance=None):
    """
    Correct and collapse the umi_per_lineage file using centroid information.

    Parameters
    ----------
    centroid_file : str
        Path to the output of `starcode_lineage_clustering_and_centroid_identification.py`.
    umi_file : str
        Path to the umi_per_lineage file.
    output_file : str
        Path to save the corrected and collapsed umi_per_lineage file.
    max_hamming_distance : int, optional
        Maximum allowable Hamming distance to filter centroids. Default is None.
    """
    logger.info("Loading centroid and umi files...")
    logger.info(centroid_file)
    logger.info(umi_file)
    centroid_df = pd.read_csv(centroid_file, sep='\t')
    umi_df = pd.read_csv(umi_file, sep='\t')

    # Optionally filter centroids by Hamming distance
    if max_hamming_distance is not None:
        logger.info(f"Filtering centroids with Hamming distance <= {max_hamming_distance}...")
        centroid_df = centroid_df[centroid_df['hamming_distance'] <= max_hamming_distance]

    max_hamming_in_df = centroid_df['hamming_distance'].max()
    logger.info(f'maximum hamming distance in centroid dataframe is {max_hamming_in_df}')

    # Step 1: Join cellbc and lineage_barcode in umi_df
    umi_df['cellbc_lineage'] = umi_df['cellbc'] + "_" + umi_df['lineage_barcode']

    # Step 2: Replace cellbc, lineage_barcode, and centroid in centroid_df
    centroid_df['cellbc_lineage'] = centroid_df['cellbc'] + "_" + centroid_df['lineage_barcode']
    centroid_df['cellbc_centroid'] = centroid_df['cellbc'] + "_" + centroid_df['centroid']

    # Step 3: Create correction dictionary
    logger.info("Creating correction dictionary from centroid dataframe...")
    correction_dict = dict(zip(centroid_df['cellbc_lineage'], centroid_df['cellbc_centroid']))

    # Step 4: Correct lineage_barcode in umi_df using correction_dict
    logger.info("Correcting lineage barcodes in umi_per_lineage file...")
    umi_df['corrected_cellbc_lineage'] = umi_df['cellbc_lineage'].map(correction_dict).fillna(umi_df['cellbc_lineage'])

    corrected_rows = umi_df[umi_df['cellbc_lineage'] != umi_df['corrected_cellbc_lineage']]
    logger.info(f'number of rows corrected {corrected_rows.shape[0]} umi_df')

    # Step 5: Separate cellbc and lineage_barcode after correction
    umi_df[['cellbc', 'lineage_barcode']] = umi_df['corrected_cellbc_lineage'].str.split("_", expand=True)

    # Collapse umi_counts based on cellbc and corrected lineage_barcode
    logger.info("Collapsing umi counts based on corrected lineage barcodes...")
    collapsed_df = (
        umi_df.groupby(['cellbc', 'lineage_barcode'], as_index=False)['umi_counts']
        .sum()
    )

    # Save the collapsed file
    collapsed_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved the collapsed file to {output_file}...")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct and collapse umi_per_lineage file using centroids.")
    parser.add_argument("--centroid_file", required=True, help="Path to the output of `starcode_lineage_clustering_and_centroid_identification.py`.")
    parser.add_argument("--umi_file", required=True, help="Path to the umi_per_lineage file.")
    parser.add_argument("--output_file", required=True, help="Path to save the corrected and collapsed umi_per_lineage file.")
    parser.add_argument("--max_hamming_distance", type=int, help="Maximum allowable Hamming distance to filter centroids. Default is None.")
    parser.add_argument("--log_folder", default=None, help="Folder to save the log file.")
    parser.add_argument("--log_config", default=None, help="Path to the logging configuration file.")
    args = parser.parse_args()

    script_name = "correct_lineages_with_starcode_and_collapse_umi"  # no file extension .py
    logger = setup_logging(script_name, args.log_config, args.log_folder)

    correct_and_collapse_umi_per_lineage(
        centroid_file=args.centroid_file,
        umi_file=args.umi_file,
        output_file=args.output_file,
        max_hamming_distance=args.max_hamming_distance
    )

