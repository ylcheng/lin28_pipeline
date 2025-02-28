#!/usr/bin/env python3

import argparse


def parse_fastq(file_path, barcode_length):
    """
    Parse a FASTQ file to extract barcodes and their associated quality scores.

    Parameters:
        file_path (str): Path to the FASTQ file.
        barcode_length (int): Length of the barcode (number of bases from the start of each read).

    Yields:
        tuple: A tuple containing:
            - str: Barcode sequence (first `barcode_length` bases of the read).
            - list: List of quality scores (Phred values) for the barcode.
    """
    if barcode_length <= 0:
        raise ValueError("Barcode length must be a positive integer.")

    with open(file_path, 'r') as fastq:
        while True:
            # Read four lines at a time (one FASTQ entry)
            identifier = fastq.readline().strip()  # Line 1: Identifier
            sequence = fastq.readline().strip()    # Line 2: Sequence
            plus = fastq.readline().strip()        # Line 3: Separator (+)
            quality = fastq.readline().strip()     # Line 4: Quality scores

            if not identifier or not sequence or not plus or not quality:
                # End of file or incomplete entry
                break

            if len(sequence) < barcode_length or len(quality) < barcode_length:
                raise ValueError(
                    f"Entry has insufficient length for barcode extraction: {identifier}"
                )

            # Extract barcode and quality scores
            barcode = sequence[:barcode_length]
            quality_scores = [ord(char) - 33 for char in quality[:barcode_length]]

            yield barcode, quality_scores


def aggregate_barcodes(parsed_fastq, barcodes=None):
    """
    Aggregate barcodes and their quality scores from a parsed FASTQ generator.

    Parameters:
        parsed_fastq (generator): Generator yielding (barcode, quality_scores).
        barcodes (dict, optional): An existing dictionary to append to.

    Returns:
        dict: A dictionary where keys are barcodes and values are lists of quality scores.
    """
    if barcodes is None:
        barcodes = {}

    for barcode, quality_scores in parsed_fastq:
        if barcode in barcodes:
            barcodes[barcode].append(quality_scores)
        else:
            barcodes[barcode] = [quality_scores]
    return barcodes


def average_quality_scores(barcodes):
    """
    Compute average quality scores for each barcode.

    Parameters:
        barcodes (dict): A dictionary where keys are barcodes and values are lists of base quality scores.

    Returns:
        dict: A dictionary where keys are barcodes and values are average quality scores per base.
    """
    averaged_barcodes = {}
    for barcode, scores in barcodes.items():
        # Compute average quality score for each base position
        average_scores = [
            sum(position_scores) / len(position_scores)
            for position_scores in zip(*scores)
        ]
        averaged_barcodes[barcode] = average_scores
    return averaged_barcodes


def save_quality_scores(barcodes, output_file):
    """
    Save barcodes and their average quality scores to a file.

    Parameters:
        barcodes (dict): A dictionary where keys are barcodes and values are average quality scores per base.
        output_file (str): Path to the output file.
    """
    with open(output_file, 'w') as f:
        for barcode, scores in barcodes.items():
            scores_str = " ".join(f"{score:.2f}" for score in scores)
            f.write(f"{barcode}\t{scores_str}\n")


def main():
    """
    Extracts barcodes from multiple FASTQ files and computes average quality scores for each barcode.

    Input:
        - FASTQ files: Multiple files containing sequencing reads in FASTQ format.

    Output:
        - A tab-separated file with barcodes and their average quality scores.
        - Example format:
            ATCG    30.00 28.50 29.50 31.00
    """
    parser = argparse.ArgumentParser(description="Extract barcodes and compute average quality scores from multiple FASTQ files.")
    parser.add_argument("fastq", nargs='+', help="Paths to the FASTQ files (space-separated for multiple files).")
    parser.add_argument("output", help="Path to the output file for quality scores.")
    parser.add_argument("--barcode-length", type=int, required=True, help="Length of the barcode (number of bases).")
    args = parser.parse_args()

    # Initialize an empty dictionary for barcodes
    barcodes = {}

    # Process each FASTQ file
    for fastq_file in args.fastq:
        print(f"Parsing FASTQ file: {fastq_file}")
        parsed_fastq = parse_fastq(fastq_file, args.barcode_length)
        barcodes = aggregate_barcodes(parsed_fastq, barcodes)

    # Compute average quality scores for each barcode
    print("Computing average quality scores...")
    averaged_barcodes = average_quality_scores(barcodes)

    # Save the results to a file
    print(f"Saving results to {args.output}")
    save_quality_scores(averaged_barcodes, args.output)

    print(f"Quality scores saved to {args.output}")


if __name__ == "__main__":
    main()

