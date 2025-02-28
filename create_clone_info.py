import pandas as pd
import argparse
import os
import logging
import logging.config
import json
from datetime import datetime

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
    print(log_file)
    with open(log_config, 'r') as config_file:
        config = json.load(config_file)
    config['handlers']['file']['filename'] = log_file
    logging.config.dictConfig(config)
    return logging.getLogger()

def create_clone_info(lineage_df):
    """
    Create a clone info DataFrame based on lineages_count_per_cellbc.

    Parameters:
    ----------
    lineage_df : pd.DataFrame
        DataFrame containing columns ['cellbc', 'lineages_count', 'lineage_barcodes'].
        lineages_count = number of unique lineage barcodes in the cellbc
        lineage_barcodes = the list of the unique lineage barcodes separated by comma

    Returns:
    -------
    pd.DataFrame
        DataFrame with columns ['cloneid', 'lineages_count', 'lineage_barcodes', 'number_of_cells'].
    """
    # Normalize lineage_barcodes for comparison
    lineage_df['sorted_lineage_barcodes'] = (
        lineage_df['lineage_barcodes']
        .apply(lambda x: ",".join(sorted(x.split(','))))  # Split, sort, and rejoin
    )

    # Aggregate data by sorted lineage_barcodes
    clone_df = (
        lineage_df.groupby('sorted_lineage_barcodes')
        .agg(
            lineages_count=('lineages_count', 'first'),  # All rows in the same group have the same count
            number_of_cells=('cellbc', 'count')  # Count the number of cells in each clone
        )
        .reset_index()
    )

    # Sort clones by number_of_cells in descending order and assign cloneid
    clone_df = clone_df.sort_values('number_of_cells', ascending=False)
    clone_df['cloneid'] = range(1, len(clone_df) + 1)

    return clone_df


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="Create clone info from lineages_count_per_cellbc file.")
    parser.add_argument("--input_file", required=True, help="Path to the input file containing lineages_count_per_cellbc data.")
    parser.add_argument("--output_file", required=True, help="Path to save the output clone info file.")
    parser.add_argument("--log_folder", default = None, help="Folder to save the log file.")
    parser.add_argument("--log_config", default = None, help="Path to the logging configuration file.")
    args = parser.parse_args()

    script_name = "create_clone_info" # no file extension .py
    logger = setup_logging( script_name, args.log_config, args.log_folder,)


    # Load input file
    logger.info(f"Loading input file: {args.input_file}")
    lineage_df = pd.read_csv(args.input_file, sep="\t")

    # Create clone info
    logger.info("Creating clone info...")
    clone_df = create_clone_info(lineage_df)

    # Save the output file
    clone_df.to_csv(args.output_file, sep="\t", index=False)
    logger.info(f"Saved clone info to: {args.output_file}")

    if 'cloneid' in clone_df.columns:
        # Map the cloneid back to the original dataframe
        logger.info(f'Map the cloneid back to lineage_df')
        lineage_df = lineage_df.merge(
            clone_df[['sorted_lineage_barcodes', 'cloneid']],
            on='sorted_lineage_barcodes',
            how='left'
        )
        lineage_df.to_csv(args.input_file, sep = "\t", index = False)
        logger.info(f'Added cloneid to {args.input_file}')


