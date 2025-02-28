from Bio import SeqIO
from collections import Counter
from datetime import datetime
import argparse
import gzip
import logging.config
import json
import os
# from lin28_pipeline.extract_lineage_barcode_data import find_flanking_sequences
from extract_lineage_barcode_data import find_flanking_sequences


def setup_logging(log_folder_path, config_folder_path, today_date):
    """
    Set up logging configuration with a dynamic log file.

    Parameters
    ----------
    log_folder_path : str
        Path to the folder where the log file will be saved.
    config_folder_path : str
        Path to the folder containing the logging configuration file.
    today_date : str
        Today's date to include in the log file name.
    """
    log_file = f"{log_folder_path}/count_cases.{today_date}.log"  # Specific log file for this script
    config_file = f"{config_folder_path}/logging_config.json"  # Path to the config file

    # Load and configure logging
    with open(config_file, "r") as config_file:
        config = json.load(config_file)
    config["handlers"]["file"]["filename"] = log_file  # Set the log file dynamically
    config["handlers"]["file"]["mode"] = "w"
    logging.config.dictConfig(config)
    return logging.getLogger()


def count_cases(r2_file, start_flanking, stop_flanking, logger):
    """
    Count the occurrences of each case (1 to 4) in the R2 FASTQ file.

    This function processes the Read 2 FASTQ file, identifying the occurrences
    of various cases based on the positions of start and stop flanking sequences.

    Parameters
    ----------
    r2_file : str
        Path to the gzipped FASTQ file for Read 2.
    start_flanking : str
        The start flanking sequence to search for.
    stop_flanking : str
        The stop flanking sequence to search for.
    logger : logging.Logger
        Logger instance for logging results.
    """
    case_counts = Counter()  # Initialize a counter for each case

    # Read the FASTQ file using SeqIO
    with gzip.open(r2_file, 'rt') as r2_stream:
        for record in SeqIO.parse(r2_stream, "fastq"):
            r2_read = str(record.seq)  # Extract the sequence from the FASTQ record

            # Find start and stop flanking sequence positions
            start_positions, stop_positions = find_flanking_sequences(r2_read, start_flanking, stop_flanking)

            # Case 1: Multiple matches for start or stop flanking sequences
            if len(start_positions) > 1 or len(stop_positions) > 1:
                case_counts["Case 1: Multiple Matches"] += 1
                continue

            # Case 2: Single match for both start and stop flanking sequences
            elif len(start_positions) == 1 and len(stop_positions) == 1:
                start_pos, stop_pos = start_positions[0], stop_positions[0]
                # if stop_pos - start_pos == len(start_flanking) + LINEAGE_BARCODE_LENGTH:
                if stop_pos - start_pos == len(start_flanking) + 20:
                    case_counts["Case 2: Single Start and Stop Match"] += 1
                else:
                    case_counts["Case 2: Invalid Region Length"] += 1
                continue

            # Case 3: Single match for the start flanking sequence only
            elif len(start_positions) == 1:
                case_counts["Case 3: Single Start Match"] += 1
                continue

            # Case 4: Single match for the stop flanking sequence only
            elif len(stop_positions) == 1:
                case_counts["Case 4: Single Stop Match"] += 1
                continue

            # Handle cases where no matches are found
            else:
                case_counts["No Match Found"] += 1

    # Log and print results
    logger.info("Case Counts:")
    for case, count in case_counts.items():
        logger.info(f"{case}: {count}")
        print(f"{case}: {count}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count cases in R2 FASTQ file based on start and stop flanking sequences.")
    parser.add_argument("--r2_file", required=True, help="Path to the gzipped FASTQ file for Read 2.")
    parser.add_argument("--start_flanking", required=True, help="Start flanking sequence.")
    parser.add_argument("--stop_flanking", required=True, help="Stop flanking sequence.")
    parser.add_argument("--log_folder_path", required=True, help="Path to the folder for logs.")
    parser.add_argument("--config_folder_path", required=True, help="Path to the folder for logging configuration.")
    args = parser.parse_args()

    today_date = datetime.now().strftime("%Y-%m-%d")  # Get today's date
    logger = setup_logging(args.log_folder_path, args.config_folder_path, today_date)  # Set up logging

    start_time = datetime.now()  # Start timing
    logger.info(f"EXECUTING count_flanking_cases.py")
    logger.info(f"--r2_file {args.r2_file}")
    logger.info(f"--start_flanking {args.start_flanking}")
    logger.info(f"--stop_flanking {args.stop_flanking}")

    count_cases(args.r2_file, args.start_flanking, args.stop_flanking, logger)  # Count cases

    end_time = datetime.now()  # End timing
    elapsed_time = end_time - start_time
    logger.info(f"COMPLETED count_flanking_cases.py; DURATION: {elapsed_time}")

