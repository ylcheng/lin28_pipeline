import argparse
import os
import gzip
import numpy as np
from Bio import SeqIO
import re
import logging.config
import json



def write_line(file, row, separator = "\t"):
    file.write(separator.join(row) + "\n")


def compute_number_of_lineage_barcode_mismatches(lineage_barcode, lineage_barcode_length = 20):
    """
    Calculate mismatches in a lineage barcode based on predefined structural rules: NSNWSNWSNWSNWSNWSNW

    This function checks specific positions in the lineage barcode:
    - Positions 1, 5, 9, 13, and 17 (0-based) must be Strong bases 'G' or 'C'.
    - Positions 3, 7, 11, 15, and 19 (0-based) must be Weak bases 'A' or 'T'.

    Any deviation at these positions is counted as a mismatch. The function returns the total
    number of mismatches.

    Parameters
    ----------
    lineage_barcode : str
        The DNA sequence representing the lineage barcode. Must be 20 bases long.

    Returns
    -------
    int
        The total number of mismatches in the barcode.
    """
    if len(lineage_barcode) == lineage_barcode_length:
        S = {'G', 'C'}
        W = {'A', 'T'}
        s_errors = sum(1 for i in [1, 5, 9, 13, 17] if lineage_barcode[i] not in S)
        w_errors = sum(1 for i in [3, 7, 11, 15, 19] if lineage_barcode[i] not in W)
        return s_errors + w_errors
    else:
        return np.nan

def compute_hamming_distance(seq1, seq2):
    """
    Calculate the Hamming distance between two sequences.

    The Hamming distance is defined as the number of positions at which
    the corresponding characters in two sequences differ. The sequences
    must be of equal length to compute the distance.

    Parameters
    ----------
    seq1 : str
        The first sequence for comparison.
    seq2 : str
        The second sequence for comparison.

    Returns
    -------
    int
        The Hamming distance between the two sequences.
    """
    if len(seq1) != len(seq2):
        return np.nan
    else:
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


# def fastq_paired_reads_generator(r1_file, r2_file, logger):
def fastq_paired_reads_generator(r1_file, r2_file):
    """
    Generate paired sequences from gzipped FASTQ files.

    This function reads two gzipped FASTQ files (Read 1 and Read 2) simultaneously,
    ensuring that paired sequences from both files are processed together. It creates
    a generator object that reads one record pair at a time when iterated over.

    Parameters
    ----------
    r1_file : str
        Path to the gzipped FASTQ file for Read 1.
    r2_file : str
        Path to the gzipped FASTQ file for Read 2.

    Yields
    ------
    tuple of str
        A tuple containing:
        - The sequence from a record in Read 1.
        - The sequence from the corresponding record in Read 2.
    """
    with gzip.open(r1_file, 'rt') as r1_stream, gzip.open(r2_file, 'rt') as r2_stream:
        r1_records = SeqIO.parse(r1_stream, 'fastq')
        r2_records = SeqIO.parse(r2_stream, 'fastq')
        for r1_record, r2_record in zip(r1_records, r2_records):
            r1_id = r1_record.id
            r2_id = r2_record.id
            if r1_id != r2_id:
                error_message = f"Record ID mismatch: R1 ID = '{r1_id}', R2 ID = '{r2_id}'"
                logger.error(error_message)
                raise ValueError(error_message)
            yield r2_id, str(r1_record.seq), str(r2_record.seq)


def find_flanking_sequences(sequence, start_flanking, stop_flanking):
    """
    Find all matches for start and stop flanking sequences in a given sequence.

    This function uses regular expressions to efficiently find the starting positions of
    all occurrences of the start and stop flanking sequences within the input sequence.

    Parameters
    ----------
    sequence : str
        The sequence in which to search for the flanking sequences.
    start_flanking : str
        The flanking sequence marking the start of the region of interest.
    stop_flanking : str
        The flanking sequence marking the end of the region of interest.

    Returns
    -------
    tuple of list of int
        A tuple containing two lists:
        - `start_positions`: The starting indices of all matches for the start flanking sequence.
        - `stop_positions`: The starting indices of all matches for the stop flanking sequence.
    """
    # Find all start flanking positions
    start_positions = [match.start() for match in re.finditer(re.escape(start_flanking), sequence)]
    # Find all stop flanking positions
    stop_positions = [match.start() for match in re.finditer(re.escape(stop_flanking), sequence)]

    return start_positions, stop_positions


def extract_sequences_based_on_start_flanking_position(start_flanking_start_pos, r2_read, start_flanking_length , stop_flanking_length , lineage_barcode_length):
    """
    Sequence structure: start_flanking-lineage_barcode-stop_flanking
    Extract flanking sequence and lineage barcode sequence from r2_read based on the known start position of the start_flanking
    Use in Case 3 and Case 1
    """
    start_flanking_end_pos = start_flanking_start_pos + start_flanking_length
    lineage_barcode_start_pos = start_flanking_end_pos
    lineage_barcode_end_pos = lineage_barcode_start_pos + lineage_barcode_length
    stop_flanking_start_pos = lineage_barcode_end_pos
    stop_flanking_end_pos = stop_flanking_start_pos + stop_flanking_length
    start_flanking = r2_read[start_flanking_start_pos:start_flanking_end_pos]
    lineage_barcode = r2_read[lineage_barcode_start_pos:lineage_barcode_end_pos]
    stop_flanking = r2_read[stop_flanking_start_pos:stop_flanking_end_pos]
    return start_flanking, lineage_barcode, stop_flanking


def extract_sequences_based_on_stop_flanking_position( stop_flanking_start_pos, r2_read, start_flanking_length , stop_flanking_length , lineage_barcode_length):
    """
    Sequence structure: start_flanking-lineage_barcode-stop_flanking
    Extract flanking sequence and lineage barcode sequence from r2_read based on the known start position of the stop_flanking
    Use in Case 4 and Case 1
    """
    stop_flanking_end_pos = stop_flanking_start_pos + stop_flanking_length
    lineage_barcode_end_pos = stop_flanking_start_pos
    lineage_barcode_start_pos = lineage_barcode_end_pos - lineage_barcode_length
    start_flanking_end_pos = lineage_barcode_start_pos
    start_flanking_start_pos = start_flanking_end_pos - start_flanking_length
    stop_flanking = r2_read[stop_flanking_start_pos:stop_flanking_end_pos]
    lineage_barcode = r2_read[lineage_barcode_start_pos:lineage_barcode_end_pos]
    start_flanking = r2_read[start_flanking_start_pos:start_flanking_end_pos]
    return start_flanking, lineage_barcode, stop_flanking

def extract_sequences_based_on_both_flanking_positions( start_flanking_start_pos, stop_flanking_start_pos, r2_read, start_flanking_length , stop_flanking_length , lineage_barcode_length):
    """
    Sequence structure: start_flanking-lineage_barcode-stop_flanking
    Extract flanking sequence and lineage barcode sequence from r2_read based on the known start positions of both start_flanking and stop_flanking
    Use in Case 2
    """
    start_flanking_end_pos = start_flanking_start_pos + start_flanking_length
    stop_flanking_end_pos = stop_flanking_start_pos + stop_flanking_length
    lineage_barcode_start_pos = start_flanking_end_pos
    lineage_barcode_end_pos = stop_flanking_start_pos
    start_flanking = r2_read[start_flanking_start_pos:start_flanking_end_pos]
    stop_flanking = r2_read[stop_flanking_start_pos:stop_flanking_end_pos]
    lineage_barcode = r2_read[lineage_barcode_start_pos:lineage_barcode_end_pos]
    return start_flanking, lineage_barcode, stop_flanking




def process_reads(r1_file, r2_file, start_flanking, stop_flanking, output_file, troubleshoot_file, lineage_barcode_length , cellbc_length , umi_length):
    """
    Process paired reads from gzipped FASTQ files and extract information for valid and invalid cases.

    This function processes paired reads from Read 1 (R1) and Read 2 (R2) FASTQ files to identify
    flanking sequences and lineage barcodes. It extract the sequence data and writes
    the results to two output files:
    - A primary output file for valid cases.
    - A troubleshooting file for ambiguous cases.

    The function matches start and stop flanking sequences in R2 and associates them with
    corresponding cell barcodes (cellbc) and unique molecular identifiers (umi) from R1.

    Expected region of interest in R2: start_flanking-lineage_barcode-stop_flanking
    Expected sequence in R1: cellbc-umi

    Parameters
    ----------
    r1_file : str
        Path to the gzipped FASTQ file containing Read 1 sequences.
    r2_file : str
        Path to the gzipped FASTQ file containing Read 2 sequences.
    start_flanking : str
        The flanking sequence marking the start of the region of interest in R2.
    stop_flanking : str
        The flanking sequence marking the end of the region of interest in R2.
    output_file : str
        Path to the output file where valid read information will be written.
    troubleshoot_file : str
        Path to the output file where invalid or ambiguous read information will be written.

    Returns
    -------
    None

    Outputs
    -------
    - output_file:
        Contains columns:
        - cellbc: The cell barcode extracted from R1.
        - umi: The unique molecular identifier extracted from R1.
        - lineage_barcode: The sequence between the start and stop flanking sequences in R2, with valid length
        - number_of_mismatch_lineage_barcode: Number of mismatches in the lineage barcode based on validation rules.
        - start_flanking: The start flanking sequence.
        - number_of_mismatch_start_flanking: Hamming distance from the expected start flanking sequence.
        - stop_flanking: The stop flanking sequence.
        - number_of_mismatch_stop_flanking: Hamming distance from the expected stop flanking sequence.

    - troubleshoot_file:
        Contains the same columns but records reads with:
        - Multiple matches for start or stop flanking sequences.
        - Match of start and stop flanking sequences, but Lineage barcodes that do not meet the length criterium.
    """
    start_flanking_length = len(start_flanking)
    stop_flanking_length = len(stop_flanking)
    # output files
    with open(output_file, 'w') as out_f, open(troubleshoot_file, 'w') as trouble_f:
        # Write headers
                # Define column headers for the output files
        headers = [
                "read_id",
            "cellbc", "umi", "lineage_barcode", "number_of_mismatch_lineage_barcode",
            "start_flanking", "number_of_mismatch_start_flanking",
            "stop_flanking", "number_of_mismatch_stop_flanking"
        ]
        # Write the headers to both output and troubleshoot files
        write_line(out_f, headers)
        write_line(trouble_f, headers)

        # Iterate over paired reads from the two FASTQ files
        for read_id, r1_read, r2_read in fastq_paired_reads_generator(r1_file, r2_file): 
            # Extract the cell barcode (first 16 bases) and UMI (next 12 bases) from Read 1
            cellbc = r1_read[:cellbc_length]
            umi = r1_read[cellbc_length: cellbc_length+umi_length]
            # Find start positions of start and stop flanking sequences in Read 2
            start_positions, stop_positions  = find_flanking_sequences(r2_read, start_flanking, stop_flanking)
            # Case 1: Handle multiple matches for start or stop flanking sequences
            if len(start_positions) > 1 or len(stop_positions) > 1:
                # Process each match for the start flanking sequence
                for start_pos in start_positions:
                    extracted_start_flanking, extracted_lineage_barcode, extracted_stop_flanking = extract_sequences_based_on_start_flanking_position( start_flanking_start_pos, r2_read, start_flanking_length , stop_flanking_length , lineage_barcode_length) 
                    number_of_mismatch_start_flanking = 0
                    number_of_mismatch_stop_flanking = compute_hamming_distance(stop_flanking, extracted_stop_flanking)
                    number_of_mismatch_lineage_barcode = compute_number_of_lineage_barcode_mismatches(extracted_lineage_barcode)
                    row = [read_id, cellbc, umi, extracted_lineage_barcode, str(number_of_mismatch_lineage_barcode), extracted_start_flanking, str(number_of_mismatch_start_flanking), extracted_stop_flanking, str(number_of_mismatch_stop_flanking)]
                    # Write the row to the troubleshooting file
                    write_line(trouble_f, row)
                # Process each match for the stop flanking sequence
                for stop_pos in stop_positions:
                    extracted_start_flanking, extracted_lineage_barcode, extracted_stop_flanking = extract_sequences_based_on_stop_flanking_position( stop_flanking_start_pos, r2_read, start_flanking_length , stop_flanking_length , lineage_barcode_length)
                    number_of_mismatch_start_flanking = compute_hamming_distance(start_flanking, extracted_start_flanking)
                    number_of_mismatch_stop_flanking = 0
                    number_of_mismatch_lineage_barcode = compute_number_of_lineage_barcode_mismatches(extracted_lineage_barcode)
                    row = [read_id, cellbc, umi, extracted_lineage_barcode, str(number_of_mismatch_lineage_barcode), extracted_start_flanking, str(number_of_mismatch_start_flanking), extracted_stop_flanking, str(number_of_mismatch_stop_flanking)]
                    # Write the row to the troubleshooting file
                    write_line(trouble_f, row)

            # Case 2: Handle a single match for both start and stop flanking sequences
            elif len(start_positions) == 1 and len(stop_positions) == 1:
                # Extract the positions of the single matches
                start_flanking_start_pos, stop_flanking_start_pos = start_positions[0], stop_positions[0]
                extracted_start_flanking, extracted_lineage_barcode, extracted_stop_flanking = extract_sequences_based_on_both_flanking_positions( start_flanking_start_pos, stop_flanking_start_pos, r2_read, start_flanking_length , stop_flanking_length , lineage_barcode_length)
                number_of_mismatch_start_flanking = 0
                number_of_mismatch_stop_flanking = 0
                if len(extracted_lineage_barcode) == lineage_barcode_length:
                    number_of_mismatch_lineage_barcode = compute_number_of_lineage_barcode_mismatches(extracted_lineage_barcode)
                    row = [read_id, cellbc, umi, extracted_lineage_barcode, str(number_of_mismatch_lineage_barcode), extracted_start_flanking, str(number_of_mismatch_start_flanking), extracted_stop_flanking, str(number_of_mismatch_stop_flanking)]
                    # Write the row to the output file
                    write_line(out_f, row)
                else:
                    # lineage_barcode length between the start_flanking and stop_flanking does not match the expected
                    number_of_mismatch_lineage_barcode = np.nan
                    row = [read_id, cellbc, umi, extracted_lineage_barcode, str(number_of_mismatch_lineage_barcode), extracted_start_flanking, str(number_of_mismatch_start_flanking), extracted_stop_flanking, str(number_of_mismatch_stop_flanking)]
                    # Write the row to the troubleshooting file
                    write_line(trouble_f, row)

            # Case 3: Handle a single match for the start flanking sequence only
            elif len(start_positions) == 1 and len(stop_positions) == 0:
                start_flanking_start_pos = start_positions[0]
                extracted_start_flanking, extracted_lineage_barcode, extracted_stop_flanking = extract_sequences_based_on_start_flanking_position( start_flanking_start_pos, r2_read, start_flanking_length , stop_flanking_length , lineage_barcode_length) 
                number_of_mismatch_start_flanking = 0
                number_of_mismatch_stop_flanking = compute_hamming_distance(stop_flanking, extracted_stop_flanking)
                number_of_mismatch_lineage_barcode = compute_number_of_lineage_barcode_mismatches(extracted_lineage_barcode)
                row = [read_id, cellbc, umi, extracted_lineage_barcode, str(number_of_mismatch_lineage_barcode), extracted_start_flanking, str(number_of_mismatch_start_flanking), extracted_stop_flanking, str(number_of_mismatch_stop_flanking)]
                # Write the row to the output file
                write_line(out_f, row)

            # Case 4: Handle a single match for the stop flanking sequence only
            elif len(start_positions) == 0 and len(stop_positions) == 1:
                stop_flanking_start_pos = stop_positions[0]
                extracted_start_flanking, extracted_lineage_barcode, extracted_stop_flanking = extract_sequences_based_on_stop_flanking_position( stop_flanking_start_pos, r2_read, start_flanking_length , stop_flanking_length , lineage_barcode_length)
                number_of_mismatch_start_flanking = compute_hamming_distance(start_flanking, extracted_start_flanking)
                number_of_mismatch_stop_flanking = 0
                number_of_mismatch_lineage_barcode = compute_number_of_lineage_barcode_mismatches(extracted_lineage_barcode)
                row = [read_id, cellbc, umi, extracted_lineage_barcode, str(number_of_mismatch_lineage_barcode), extracted_start_flanking, str(number_of_mismatch_start_flanking), extracted_stop_flanking, str(number_of_mismatch_stop_flanking)]
                # Write the row to the output file
                write_line(out_f, row)





def parse_arguments():
    parser = argparse.ArgumentParser(description="Process paired reads from gzipped FASTQ files to extract cellbc, umi, lineage_barcode, start_flanking, and stop_flanking. Compute the mismatched between the extracted sequences and expected sequences.")

    parser.add_argument("--r1_file", required=True, help="Path to the gzipped FASTQ file for Read 1.")
    parser.add_argument("--r2_file", required=True, help="Path to the gzipped FASTQ file for Read 2.")
    parser.add_argument("--start_flanking", default='GAAACACCG', help="Start flanking sequence to search for in Read 2.")
    parser.add_argument("--stop_flanking", default='GTTTTAGAG', help="Stop flanking sequence to search for in Read 2.")
    parser.add_argument("--output_file", required=True, help="Path to the output file for valid reads.")
    parser.add_argument("--troubleshoot_file", required=True, help="Path to the output file for troubleshooting ambiguous reads.")
    parser.add_argument("--lineage_barcode_length", type=int, default=20, help="Length of the lineage barcode (default: 20).")
    parser.add_argument("--cellbc_length", type=int, default=16, help="Length of the cell barcode in Read 1 (default: 16).")
    parser.add_argument("--umi_length", type=int, default=12, help="Length of the UMI in Read 1 (default: 12).")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    # Setup logging
    log_config = os.getenv("LOG_CONFIG")
    log_file = os.getenv("LOG_FILE")
    with open(log_config, "r") as config_file:
        config = json.load(config_file)
    config["handlers"]["file"]["filename"] = log_file
    logging.config.dictConfig(config)
    logger = logging.getLogger()

    process_reads(
        r1_file=args.r1_file,
        r2_file=args.r2_file,
        start_flanking=args.start_flanking,
        stop_flanking=args.stop_flanking,
        output_file=args.output_file,
        troubleshoot_file=args.troubleshoot_file,
        lineage_barcode_length=args.lineage_barcode_length,
        cellbc_length=args.cellbc_length,
        umi_length=args.umi_length,
            )





