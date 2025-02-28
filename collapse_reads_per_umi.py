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


def collapse_reads(input_file, output_file, separator = '\t'):
    """
    Collapse rows with the same 'cellbc', 'umi', and 'lineage_barcode' into unique molecules.

    Parameters
    ----------
    input_file : str
        Path to the input TSV file containing the filtered reads.
    output_file : str
        Path to the output TSV file where collapsed data will be saved.
    separator : str, optional
        Separator used in the input file, default is tab ("\t").
    """
    logger.info(f"Loading data {input_file}...")
    df = pd.read_csv(input_file, sep = separator)

    # Verify requried columns exist
    required_columns = ['cellbc', 'umi', 'lineage_barcode']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        logger.error(f'Missing required columns: {", ".join(missing_columns)}')
        raise ValueError(f'Required columns not found: {", ".join(missing_columns)}')

    # Combine columns to identify unique molecules
    logger.info(f'Concatenate cellbc, umi, lineage_barcode to mark molecule_id')
    logger.info("Collapsing reads per umi...")
    df['molecule_id'] = df[['cellbc', 'umi','lineage_barcode']].agg('_'.join, axis = 1)
    # Group by 'molecule_id' and count the number of reads 
    collapsed_df = (df.groupby('molecule_id').size().reset_index(name = 'reads_count'))
    # Split 'molecule_id' back into original columns
    collapsed_df[['cellbc', 'umi', 'lineage_barcode']] = collapsed_df['molecule_id'].str.split('_', expand = True)
    # Reorder columns and drop 'molecule_id
    collapsed_df = collapsed_df[['cellbc', 'umi', 'lineage_barcode', 'reads_count']]
    
    logger.info(f"Saving collapsed data to {output_file}...")
    collapsed_df.to_csv(output_file, sep = separator, index = False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collapse rows with the same cell barcode, umi, and lineage barcode into unique molecules and count reads.")
    parser.add_argument("--input_file", required=True, help="Path to the input TSV file.")
    parser.add_argument("--output_file", required=True, help="Path to the output TSV file.")
    parser.add_argument("--separator", default="\t", help="Separator used in the input and output files (default: tab).")
    parser.add_argument("--log_folder", default = None, help="Folder to save the log file.")
    parser.add_argument("--log_config", default = None, help="Path to the logging configuration file.")
    args = parser.parse_args()

    script_name = "collapse_reads_per_umi"
    logger = setup_logging( script_name, args.log_config, args.log_folder,)

    # Collapse reads per umi
    collapse_reads(args.input_file, args.output_file, args.separator)







