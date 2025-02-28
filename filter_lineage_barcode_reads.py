import pandas as pd
import logging
import logging.config
import json
import argparse
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
    print(log_folder)
    with open(log_config, 'r') as config_file:
        config = json.load(config_file)
    config['handlers']['file']['filename'] = log_file
    logging.config.dictConfig(config)
    return logging.getLogger()


def apply_mismatch_filters(df_valid, max_lineage_mismatch, max_start_flanking_mismatch, max_stop_flanking_mismatch):
    """
    Apply mismatch filters to the lineage barcode reads.
    """
    lineage_cond = df_valid['number_of_mismatch_lineage_barcode'] <= max_lineage_mismatch
    start_cond = df_valid['number_of_mismatch_start_flanking'] <= max_start_flanking_mismatch
    stop_cond = df_valid['number_of_mismatch_stop_flanking'] <= max_stop_flanking_mismatch
    return df_valid[lineage_cond & start_cond & stop_cond]


def correct_and_filter_cell_barcodes(df_valid, hamming_df, rna_cell_barcodes, max_hamming_distance, ):
    """
    Correct and filter cell barcodes based on cell barcodes from RNA dataset and maximum allowable Hamming distance.
    """
    if max_hamming_distance == 0:
        # only keep the cell barcodes in lineage dataset if they have exact match with the cell barcodes in RNA dataset
        logger.info("Max Hamming distance is 0. Skipping correction and keeping only exact matches.")
        return df_valid[df_valid['cellbc'].isin(rna_cell_barcodes)]

    else:
        # Step 1: Filter hamming_df for max_hamming_distance
        hamming_filtered = hamming_df[hamming_df['HammingDistance'] <= max_hamming_distance]
        # filter out Query that have more than one match of Reference that are same hamming distance away
        # if hamming distances are different, kept the Query-Reference pair that has the lowest hamming distance

        # Step 2: Create a correction dictionary
        correction_dict = dict(zip(hamming_filtered['Query'], hamming_filtered['Reference']))
        # Step 3: Filter df_valid for exact matches and corrected matches
        exact_match_df = df_valid[df_valid['cellbc'].isin(rna_cell_barcodes)]
        correction_match_df = df_valid[df_valid['cellbc'].isin(correction_dict.keys())]
        correction_match_df['cellbc'] = correction_match_df['cellbc'].map(correction_dict)
        combined_df = pd.concat([exact_match_df, correction_match_df], ignore_index = True)
        return combined_df

def correct_and_filter_cell_barcodes(df_valid, hamming_df, rna_cell_barcodes, max_hamming_distance):
    """
    Correct and filter cell barcodes based on cell barcodes from RNA dataset and maximum allowable Hamming distance.

    Parameters:
    ----------
    df_valid : pd.DataFrame
        DataFrame containing valid lineage barcodes.
    hamming_df : pd.DataFrame
        DataFrame containing Hamming distances between Queries and References.
    rna_cell_barcodes : set
        Set of RNA cell barcodes (reference barcodes).
    max_hamming_distance : int
        Maximum allowable Hamming distance for correction.

    Returns:
    -------
    pd.DataFrame
        DataFrame with corrected and filtered lineage barcodes.
    """
    if max_hamming_distance == 0:
        # Only keep the cell barcodes in lineage dataset if they have an exact match with the cell barcodes in RNA dataset
        logger.info("Max Hamming distance is 0. Skipping correction and keeping only exact matches.")
        return df_valid[df_valid['cellbc'].isin(rna_cell_barcodes)]

    else:
        # Step 1: Filter hamming_df for max_hamming_distance
        hamming_filtered = hamming_df[hamming_df['HammingDistance'] <= max_hamming_distance]

        # Step 2: Resolve multiple matches for the same Query

        # filter function to resolve multiple matches for the same Query
        def filter_group(group, removed_queries):
            # Find the smallest hamming distance in the group
            min_distance = group['HammingDistance'].min()
            # Count how many rows have the smalles hamming distance
            min_distance_count = (group['HammingDistance'] == min_distance).sum()
            if min_distance_count > 1:
                # If more than one row has the smallest Hamming distance, remove all rows for this Query
                removed_queries.append(group['Query'].iloc[0]) # Log the Query being removed
                return pd.DataFrame() # Return an empty DataFrame for this group
            else:
                # Otherwise, keep only the row with the smallest Hamming distance
                return group[group['HammingDistance'] == min_distance]

        # Separate rows with unique Query from non-unique Query
        query_counts = hamming_filtered['Query'].value_counts()
        unique_queries = query_counts[query_counts == 1].index
        non_unique_queries = query_counts[query_counts > 1].index
        # Extract unique and non-unique rows
        unique_rows = hamming_filtered[hamming_filtered['Query'].isin(unique_queries)]
        non_unique_rows = hamming_filtered[hamming_filtered['Query'].isin(non_unique_queries)]
        # Apply filter_group only to non_unique_queries
        removed_queries = [] # to log removed queries
        filtered_non_unique_rows = (
                non_unique_rows
                .sort_values(by = ['Query', 'HammingDistance'])
                .groupby('Query', group_keys = False)
                .apply(lambda group: filter_group(group, removed_queries))
                )
        # Combine filtered rows
        hamming_filtered = pd.concat([unique_rows, filtered_non_unique_rows], ignore_index = True)


        logger.info(f"Removed {len(removed_queries)} queries due to ties in the smallest Hamming distance: {removed_queries}")

        # Step 3: Create a correction dictionary
        correction_dict = dict(zip(hamming_filtered['Query'], hamming_filtered['Reference']))

        # Step 4: Filter df_valid for exact matches and corrected matches
        exact_match_df = df_valid[df_valid['cellbc'].isin(rna_cell_barcodes)]
        correction_match_df = df_valid[df_valid['cellbc'].isin(correction_dict.keys())]
        correction_match_df['cellbc'] = correction_match_df['cellbc'].map(correction_dict)

        # Combine exact matches and corrected matches
        combined_df = pd.concat([exact_match_df, correction_match_df], ignore_index=True)
        return combined_df



def main(valid_file, hamming_file, rna_file, max_hamming_distance,
         max_lineage_mismatch, max_start_flanking_mismatch, max_stop_flanking_mismatch,
         output_file):
    """
    Main function to filter lineage barcode reads based on various criteria.

    Filtering steps:
    ----------------
    1. **cell barcode: Exact Match Filtering (Always Performed)**:
       - Retains rows from `valid_file` where the `cellbc` (cell barcode) matches exactly with
         cell barcodes from the RNA dataset (`rna_file`).

    2. **cell barcode: Hamming Distance Correction and Filtering** (if `max_hamming_distance > 0`):
       - Uses the precomputed Hamming distance file (`hamming_file`) to find lineage cell barcodes
         that are within the specified Hamming distance from RNA cell barcodes.
       - Constructs a correction dictionary to map lineage cell barcodes to their closest RNA
         barcode match based on Hamming distance.
       - Updates the `cellbc` field in `valid_file` with corrected RNA barcodes for matched rows.
       - Combines exact matches and corrected barcodes into a single DataFrame.

    3. **lineage barcode and flanking: Mismatch Filtering**:
       - Applies thresholds to the following columns in `valid_file`:
         - `number_of_mismatch_lineage_barcode`: Number of mismatches in the lineage barcode.
         - `number_of_mismatch_start_flanking`: Mismatches in the start flanking sequence.
         - `number_of_mismatch_stop_flanking`: Mismatches in the stop flanking sequence.
       - Retains only rows that satisfy all mismatch thresholds.

    Output:
    -------
    - Saves the filtered DataFrame to `output_file`, with the filename dynamically
      updated to include the applied thresholds (e.g., `.cellbc1start1lineage2stop1.tsv`).

    Parameters:
    -----------
    valid_file : str
        Path to the input file containing valid lineage barcode reads.
        The file must have the following columns: `cellbc`, `number_of_mismatch_lineage_barcode`,
        `number_of_mismatch_start_flanking`, `number_of_mismatch_stop_flanking`.

    hamming_file : str, optional
        Path to the precomputed Hamming distances file. Must include `Query`, `Reference`,
        and `HammingDistance` columns. Required if `max_hamming_distance > 0`.

    rna_file : str
        Path to the RNA dataset file containing valid RNA cell barcodes. The file should
        have no header, and RNA barcodes should be in the first column.

    max_hamming_distance : int
        Maximum allowable Hamming distance for cell barcode correction. If `0`, only exact
        matches with RNA barcodes are retained, and `hamming_file` is ignored.

    max_lineage_mismatch : int
        Maximum allowable mismatches for the lineage barcode.

    max_start_flanking_mismatch : int
        Maximum allowable mismatches for the start flanking sequence.

    max_stop_flanking_mismatch : int
        Maximum allowable mismatches for the stop flanking sequence.

    output_file : str
        Path to save the filtered and corrected lineage barcode data. The filename will
        include a suffix summarizing the applied thresholds.

    """



    logger.info(f"Loading valid and rna files, {valid_file}, {rna_file}...")
    df_valid = pd.read_csv(valid_file, sep = '\t')
    # logger.info(f'Number of reads in valid_file: {len(df_valid)}')
    logger.info(f'Number of reads and columns in valid_file: {df_valid.shape}')

    rna_cell_barcodes = set(pd.read_csv(rna_file, sep = '\t',  header = None).iloc[:, 0])

    if max_hamming_distance > 0:
        logger.info("Correcting and filtering cell barcodes...")
        logger.info(f'Loading hamming_file: {hamming_file}...')
        hamming_df = pd.read_csv(hamming_file, sep = '\t', compression='bz2')
        filtered_df = correct_and_filter_cell_barcodes(df_valid, hamming_df, rna_cell_barcodes, max_hamming_distance, )
    else:
        logger.info("Filtering only exactly matched cell barcodes...")
        filtered_df = df_valid[ df_valid['cellbc'].isin(rna_cell_barcodes) ]

    num_cellbc_in_rna = len(rna_cell_barcodes)
    num_cellbc_before_filter = len(df_valid['cellbc'].unique())
    num_cellbc_after_filter = len(filtered_df['cellbc'].unique())
    logger.info(f'Number of cell barcodes in RNA dataset: {num_cellbc_in_rna}')
    logger.info(f'Number of cell barcodes in lineage dataset before cell barcode filtering and correction: {num_cellbc_before_filter}')
    logger.info(f'Number of cell barcodes in lineage dataset after cell barcode filtering and correction: {num_cellbc_after_filter}')

    # Step 2: Apply mismatch thresholds
    logger.info('Applying mismatch filters on lineage barcode and flanking...')
    logger.info(f'Number of reads before mismatch filter: {filtered_df.shape[0]}')
    filtered_df = apply_mismatch_filters( filtered_df, max_lineage_mismatch, max_start_flanking_mismatch, max_stop_flanking_mismatch)
    logger.info(f'Number of reads after mismatch filter: {filtered_df.shape[0]}')
    number_cellbc_after_mismatch_filters = len(filtered_df['cellbc'].unique())
    logger.info(f'Number of cell barcodes after mismatch filter: {number_cellbc_after_mismatch_filters}')

    # Troubleshooting, suspicious that I duplicated the dataframe 
    # Check for duplicates
    num_duplicates = filtered_df.duplicated().sum()
    if num_duplicates > 0:
        logger.info(f"Number of duplicated rows: {num_duplicates}")
        logger.info(f"DataFrame before removing duplicates: {filtered_df.shape[0]}")
        logger.info(filtered_df.shape)  # Log a preview of the DataFrame before filtering

        # Remove duplicates
        filtered_df = filtered_df.drop_duplicates()

        logger.info(f"DataFrame after removing duplicates: {filtered_df.shape}")
    else:
        logger.info(f'Number of reads and columns in filtered_df: {filtered_df.shape}')



    # Step 3: Save the filtered results
    output_suffix = (f".cellbc{max_hamming_distance}"
                     f"start{max_start_flanking_mismatch}"
                     f"lineage{max_lineage_mismatch}"
                     f"stop{max_stop_flanking_mismatch}")
    output_file = output_file.replace(".tsv", f"{output_suffix}.tsv")
    filtered_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Filtered results saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and correct lineage barcodes based on RNA barcodes and mismatch thresholds.")
    parser.add_argument("--valid_file", required=True, help="Path to valid_lineage_barcode_reads.tsv.")
    parser.add_argument("--hamming_file", help="Path to precomputed cell barcodes Hamming distances file, in Bzip2 format")
    parser.add_argument("--rna_file", help="Path to RNA dataset file containing cell barcodes. Default: no header, separator is '\t', and in the first column")
    parser.add_argument("--max_hamming_distance", type=int, required=True, help="Maximum allowable Hamming distance.")
    parser.add_argument("--max_lineage_mismatch", type=int, required=True, help="Maximum allowable mismatches for lineage barcode.")
    parser.add_argument("--max_start_flanking_mismatch", type=int, required=True, help="Maximum mismatches for start flanking.")
    parser.add_argument("--max_stop_flanking_mismatch", type=int, required=True, help="Maximum mismatches for stop flanking.")
    parser.add_argument("--output_file", required=True, help="Path to save the filtered results.")
    parser.add_argument("--log_folder", default = None, help="Folder to save the log file.")
    parser.add_argument("--log_config", default = None, help="Path to the logging configuration file.")

    args = parser.parse_args()

    script_name = "filter_lineage_barcode_reads"
    logger = setup_logging( script_name, args.log_config, args.log_folder,)

    main(
        valid_file=args.valid_file,
        hamming_file=args.hamming_file,
        rna_file=args.rna_file,
        max_hamming_distance=args.max_hamming_distance,
        max_lineage_mismatch=args.max_lineage_mismatch,
        max_start_flanking_mismatch=args.max_start_flanking_mismatch,
        max_stop_flanking_mismatch=args.max_stop_flanking_mismatch,
        output_file=args.output_file,
    )











