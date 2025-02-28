import pandas as pd
import argparse
import logging
import logging.config
import json
import os

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


def collapse_molecules_and_count(input_file, output_file, separator = '\t'):
    """
    Collapse and count unique UMIs for each cellbc and lineage_barcode combination.

    Parameters
    ----------
    input_file : str
        Path to the input TSV file.
    output_file : str
        Path to the output TSV file where molecule counts will be saved.
    separator : str, optional
        Separator used in the input file, default is tab ("\t").
    """
    logger.info(f'Loading data {input_file}...')
    df = pd.read_csv(input_file, sep = separator)

    # Verify required columns exist
    required_columns = ['cellbc', 'umi', 'lineage_barcode']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        logger.error(f'Missing required columns: {", ".join(missing_columns)}')
        raise ValueError(f'Required columns not found: {", ".join(missing_columns)}')

    # Group by 'cellbc' and 'lineage_barcode' and count unique 'umi'
    logger.info(f'Group by  cellbc and lineage_barcode to count unique UMIs...')
    grouped = ( df.groupby(['cellbc', 'lineage_barcode'])['umi'].nunique().reset_index(name='umi_counts'))

    logger.info(f'Saving results {output_file}...')
    grouped.to_csv(output_file, sep=separator, index = False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collapse and count unique UMIs for each cell barcode and lineage barcode combination.")
    parser.add_argument("--input_file", required=True, help="Path to the input TSV file.")
    parser.add_argument("--output_file", required=True, help="Path to the output TSV file.")
    parser.add_argument("--separator", default="\t", help="Separator used in the input and output files (default: tab).")
    parser.add_argument("--log_folder", default = None, help="Folder to save the log file.")
    parser.add_argument("--log_config", default = None, help="Path to the logging configuration file.")
    args = parser.parse_args()

    script_name = "collapse_umi_per_lineage" # no file extension .py
    logger = setup_logging( script_name, args.log_config, args.log_folder,)

    collapse_molecules_and_count(args.input_file, args.output_file, args.separator)



