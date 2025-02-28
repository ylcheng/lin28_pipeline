# version 2
import argparse
import logging
import subprocess
from Bio import SeqIO
import gzip
import os
from glob import glob

import logging.config
import json


def find_fastq_files(folder_path, fastq_pattern, R1_pattern, R2_pattern):
    """
    Find all FASTQ files in the folder and return sorted lists of R1 and R2 files.
    """
    # Find all files matching the fastq_pattern
    fastq_files = glob(os.path.join(folder_path, fastq_pattern))

    # Check if any FASTQ files are found
    if not fastq_files:
        logger.error(f"No files matching the pattern '{fastq_pattern}' were found in the folder: {folder_path}")
        raise FileNotFoundError(f"No files matching the pattern '{fastq_pattern}' found in folder: {folder_path}")
    else: 
        logger.info(f'{len(fastq_files)} fastq.gz files found: {fastq_files}')

    r1_files = sorted([f for f in fastq_files if R1_pattern in f])
    r2_files = sorted([f for f in fastq_files if R2_pattern in f])

    r1_files_num = len(r1_files) 
    r2_files_num = len(r2_files)
    if r1_files_num != r2_files_num:
        logger.error(f"Mismatch in the number of R1 and R2 files:{r1_files_num}, {r2_files_num})")
        logger.error(f"{r1_files}")
        logger.error(f"{r2_files}")
        raise ValueError("Mismatch in the number of R1 and R2 files.")

    for i in range(len(r1_files)):
        r1_prefix = r1_files[i].split(R1_pattern)[0]
        r2_prefix = r2_files[i].split(R2_pattern)[0]
        if r1_prefix != r2_prefix:
            logger.error(f"Mismatch in the sorted R1 and R2 files: {r1_files}, {r2_files}")
            raise ValueError(f"Mismatch in the sorted R1 and R2 files: {r1_files}, {r2_files}")

    logger.info(f"Found {r1_files_num} R1 {r1_files}")
    logger.info(f"Found {r2_files_num} R2 {r2_files}")
    return r1_files, r2_files


def subset_reads_with_flanking_sequences(folder_path, fastq_pattern, R1_pattern, R2_pattern, start_sequence, stop_sequence, output_fastq_prefix):
    """
    Subset reads containing the start or stop flanking sequences in Read 2
    and save the corresponding reads in Read 1 and Read 2 to new FASTQ files.
    Process a folder of FASTQ files, subset matching reads, and directly write combined results.
    """
    # Find all FASTQ files in the folder and separate them into R1 and R2 file lists
    r1_files, r2_files = find_fastq_files(folder_path, fastq_pattern, R1_pattern, R2_pattern)

    # Generate output file paths for combined R1 and R2 files 
    combined_r1_output = os.path.join(folder_path, f'{output_fastq_prefix}_R1_flanking_detected.fastq.gz')
    combined_r2_output = os.path.join(folder_path, f'{output_fastq_prefix}_R2_flanking_detected.fastq.gz')

    # Open combined output files for writing
    with gzip.open(combined_r1_output, "wt") as combined_r1_file, gzip.open(combined_r2_output, "wt") as combined_r2_file:
        
         # Initialize counters for total reads and retained reads
         total_reads_number = 0
         retained_reads_number = 0

         # Loop through each pair of R1 and R2 files
         for read1_file, read2_file in zip(r1_files, r2_files):
             total_reads_num = 0
             retained_reads_num = 0
             logger.info(f'Subsetting: {read1_file} and {read2_file}')

             # Open the input R1 and R2 files for reading 
             with gzip.open(read1_file, "rt") as read1_input_file, gzip.open(read2_file, "rt") as read2_input_file:
                 # Parse the input files into records
                 read1_records = SeqIO.parse(read1_input_file, "fastq")
                 read2_records = SeqIO.parse(read2_input_file, "fastq")

                 # Process paired records from R1 and R2 files
                 for read1_record, read2_record in zip(read1_records, read2_records):
                     total_reads_number+=1
                     total_reads_num+=1
                     # Check if the R2 sequence contains the start or stop flanking sequence
                     if start_sequence in str(read2_record.seq) or stop_sequence in str(read2_record.seq):
                         retained_reads_number+=1
                         retained_reads_num+=1
                         # Write the matching R1 and R2 records to the combined output files
                         SeqIO.write(read1_record, combined_r1_file, "fastq")
                         SeqIO.write(read2_record, combined_r2_file, "fastq") 
             logger.info(f'total reads: {total_reads_num}, retained reads: {retained_reads_num}')
                

    # Log summary statistics
    logger.info(f"Total reads processed: {total_reads_number}")
    logger.info(f"Reads retained with flanking sequences: {retained_reads_number}")
    logger.info(f"Fraction of reads with flanking sequences: {retained_reads_number / total_reads_number:.4f}")
    logger.info(f"Combined R1 output saved to: {combined_r1_output}")
    logger.info(f"Combined R2 output saved to: {combined_r2_output}")


if __name__ == "__main__":

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Subset reads containing specific flanking sequences.")
    parser.add_argument("--folder_path", required=True, help="Path to input Read1 and Read 2 FASTQ files (gzipped).")
    parser.add_argument("--output_fastq_prefix", required=True, help="Prefix for the output combined FASTQ files retained reads with flanking sequences.")
    parser.add_argument("--start_sequence", default="GAAACACCG", help="Start flanking sequence to search for (default: GAAACACCG).")
    parser.add_argument("--stop_sequence", default="GTTTTAGAG", help="Stop flanking sequence to search for (default: GTTTTAGAG).")
    parser.add_argument("--fastq_pattern", default="*fastq.gz", help="Pattern to search for FASTQ files in the folder (default: '*fastq.gz').")
    parser.add_argument("--R1_pattern", default="_R1_", help="Pattern to identify R1 files (default: '_R1_').")
    parser.add_argument("--R2_pattern", default="_R2_", help="Pattern to identify R2 files (default: '_R2_').")
    
    args = parser.parse_args()

    # Setup logging
    # Retrieve environment variables
    log_config = os.getenv("LOG_CONFIG")
    log_file = os.getenv("LOG_FILE")
    # Load and configure logging
    with open(log_config, "r") as config_file:
        config = json.load(config_file)
    config["handlers"]["file"]["filename"] = log_file
    logging.config.dictConfig(config)
    logger = logging.getLogger()


    subset_reads_with_flanking_sequences(
            args.folder_path,
            args.fastq_pattern,
            args.R1_pattern,
            args.R2_pattern,
            args.start_sequence,
            args.stop_sequence,
            args.output_fastq_prefix,
            )
