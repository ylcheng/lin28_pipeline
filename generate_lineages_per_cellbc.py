import pandas as pd
import argparse
import logging
import logging.config
import json
import os
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

def process_umi_per_lineage(input_file, output_file):
    """
    Process umi_per_lineage.tsv to generate lineage_per_cellbc.tsv.
    
    Parameters:
    ----------
    input_file : str
        Path to the input umi_per_lineage.tsv file containing columns 'cellbc', 'lineage_barcode', 'umi_counts'.
    output_file : str
        Path to save the output lineage_per_cellbc.tsv file.
    """
    # Read the input file
    logger.info(f'Processing input {input_file}...')
    df = pd.read_csv(input_file, sep="\t")

    # Group by 'cellbc' and compute lineage counts and collect lineage barcodes
    grouped = (
        df.groupby("cellbc")["lineage_barcode"]
        .apply(lambda x: x.unique())  # Get unique lineage barcodes
        .reset_index()
    )
    
    # Calculate lineage counts and prepare the lineage_barcode list as a comma-separated string
    grouped["lineages_count"] = grouped["lineage_barcode"].apply(len)
    grouped["lineage_barcodes"] = grouped["lineage_barcode"].apply(lambda x: ",".join(x))
    grouped = grouped[["cellbc", "lineages_count", "lineage_barcodes"]]

    # Save the result to output file
    grouped.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Output saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate lineage_per_cellbc.tsv from umi_per_lineage.tsv.")
    parser.add_argument("--input_file", required=True, help="Path to the input umi_per_lineage.tsv file.")
    parser.add_argument("--output_file", required=True, help="Path to save the output lineage_per_cellbc.tsv file.")
    parser.add_argument("--log_folder", default = None, help="Folder to save the log file.")
    parser.add_argument("--log_config", default = None, help="Path to the logging configuration file.")
    args = parser.parse_args()

    script_name = "generate_lineages_per_cellbc" # no file extension .py
    logger = setup_logging( script_name, args.log_config, args.log_folder,)

    process_umi_per_lineage(args.input_file, args.output_file)

