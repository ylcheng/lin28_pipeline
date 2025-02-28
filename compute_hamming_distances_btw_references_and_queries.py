# version 2
from datetime import datetime
import pandas as pd
import numpy as np
from itertools import product
from Bio import SeqIO
from scipy.spatial import KDTree
import logging
import os
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
    Compute the Hamming distance between two sequences.
    """
    if len(seq1) != len(seq2):
        return np.nan
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def load_sequences(file_path, column_index, separator="\t", header="infer", unique=True):
    """
    Load sequences from a specific column in a file.
    """
    df = pd.read_csv(file_path, sep=separator, header=header)
    sequences = df.iloc[:, column_index].tolist()
    return set(sequences) if unique else sequences


def split_barcode(barcode, d=1):
    """
    Split the barcode into d+1 segments.
    """
    split_size = len(barcode) // (d + 1)
    return [barcode[i:i + split_size] for i in range(0, len(barcode), split_size)]


def filter_by_segment_match(reference_sequences, query_sequences, max_hamming_distance):
    """
    Filter (kept) query sequences that have at least one matching segment with reference sequences.
    if max_hamming is 3, and we split seq into 3+1 = 4 segment, for the query to potentially have max_hamming 3 with the reference, then query must have one segment matched with the refence corresponding segment. When we have no matching segment, that mean the query sequence is at least 3 hamming distance from reference sequences
    we always split query seq to d + 1 segments, so always looking for at least one matching segment with corresponding segments of reference sequences
    """
    reference_segments = [{split_barcode(seq, max_hamming_distance)[i] for seq in reference_sequences} for i in range(max_hamming_distance + 1)]

    filtered_queries = []
    for sequence in set(query_sequences):
        segments = split_barcode(sequence, max_hamming_distance)
        # looking for match of corresponding segments between query and reference
        match = any(segments[i] in reference_segments[i] for i in range(max_hamming_distance + 1))
        if match:
            filtered_queries.append(sequence)
    return set(filtered_queries)



def main(reference_file, reference_column, query_file, query_column, output_file,
         reference_separator, query_separator, reference_header, query_header,
         max_hamming_distance, segment_match):
    """
    Main function to compute Hamming distances with optional segment match filtering.
    """
    # Handle 'None' as a valid input for header
    reference_header = None if reference_header == "None" else reference_header
    query_header = None if query_header == "None" else query_header

    # Load sequences
    reference_sequences = load_sequences(reference_file, reference_column, reference_separator, reference_header)
    query_sequences = load_sequences(query_file, query_column, query_separator, query_header)
    logger.info(f"Number of unique sequences in reference set: {len(reference_sequences)}")
    logger.info(f"Number of unique sequences in query set: {len(query_sequences)}")

    # Remove exact matches from query set
    exact_matches = reference_sequences.intersection(query_sequences)
    logger.info(f"Number of exact matches removed from query set: {len(exact_matches)}")

    remaining_query_sequences = query_sequences.difference(exact_matches)

    if remaining_query_sequences:
        logger.info(f"Number of sequences remaining in query set after removing exact matches: {len(remaining_query_sequences)}")

        # Segment match filtering
        if segment_match:
            logger.info(f"Performing segment match filtering with max Hamming distance {max_hamming_distance}")
            remaining_query_sequences = filter_by_segment_match(reference_sequences, remaining_query_sequences, max_hamming_distance)
            output_file = output_file.replace(".tsv", f".max_hamming_{max_hamming_distance}.tsv")
            logger.info(f"Number of unique sequences in query set after segment match filtering: {len(query_sequences)}")
        else:
            output_file = output_file.replace(".tsv", ".all_pairwise.tsv")

        # Compute Hamming distances
        distances = []
        for ref_seq, query_seq in product(reference_sequences, remaining_query_sequences):
            hamming_distance = compute_hamming_distance(ref_seq, query_seq)
            if segment_match and hamming_distance > max_hamming_distance:
                continue
            distances.append({'Reference': ref_seq, 'Query': query_seq, 'HammingDistance': hamming_distance})

        # Save results
        df_distances = pd.DataFrame(distances)
        num_ref = len(df_distances['Reference'].unique())
        num_query = len(df_distances['Query'].unique())
        logger.info(f"Number of unique sequences in reference set in output: {num_ref}")
        logger.info(f"Number of unique sequences in query set in output: {num_query}")
        df_distances.to_csv(output_file, index=False, sep='\t', compression='bz2')
        logger.info(f"Results saved  in compressed format to {output_file}")


if __name__ == "__main__":
    import argparse
    from datetime import datetime

    parser = argparse.ArgumentParser(description="Compute pairwise Hamming distances with optional segment filtering.")
    parser.add_argument("--reference_file", required=True, help="Path to the reference file.")
    parser.add_argument("--reference_column", default = 0,  type=int, help="Column index (0-based) in the reference file.")
    parser.add_argument("--query_file", required=True, help="Path to the query file.")
    parser.add_argument("--query_column", default = 1, type=int, help="Column index (0-based) in the query file.")
    parser.add_argument("--output_file", required=True, help="Path to save the output TSV file.")
    parser.add_argument("--log_folder", default = None, help="Folder to save the log file.")
    parser.add_argument("--log_config", default = None, help="Path to the logging configuration file.")
    parser.add_argument("--reference_separator", default="\t", help="Separator for reference file. Default: '\\t'.")
    parser.add_argument("--query_separator", default="\t", help="Separator for query file. Default: '\\t'.")
    parser.add_argument("--reference_header", default=None, help="Header for reference file. Use None for no header.")
    parser.add_argument("--query_header", default="infer", help="Header for query file.")
    parser.add_argument("--max_hamming_distance", type=int, default=1, help="Maximum allowed Hamming distance for segment filtering.")
    parser.add_argument("--segment_match", action='store_true', help="Enable segment match filtering instead of all pairwise computation.")

    args = parser.parse_args()

    # Setup logging
    script_name = "compute_hamming_distances_btw_references_and_queries"
    logger = setup_logging( script_name, args.log_config, args.log_folder,)

    start_time = datetime.now()
    logger.info(f"EXECUTING {script_name}.py")
    logger.info(f"Arguments: {args}")
    main(
        reference_file=args.reference_file,
        reference_column=args.reference_column,
        query_file=args.query_file,
        query_column=args.query_column,
        output_file=args.output_file,
        reference_separator=args.reference_separator,
        query_separator=args.query_separator,
        reference_header=args.reference_header,
        query_header=args.query_header,
        max_hamming_distance=args.max_hamming_distance,
        segment_match=args.segment_match
    )
    end_time = datetime.now()
    logger.info(f"COMPLETED {script_name}.py in {end_time - start_time}")

