
import argparse


def hamming_distance(s1, s2):
    """
    Compute the Hamming distance between two strings.

    Parameters:
        s1 (str): The first string.
        s2 (str): The second string.

    Returns:
        int: The Hamming distance between the two strings.
    """
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


def error_probability(quality_score):
    """
    Convert quality score (Phred) to error probability.

    Parameters:
        quality_score (int): The quality score of a base (Phred scale).

    Returns:
        float: The probability of an error for the given quality score.
    """
    return 10 ** (-quality_score / 10)


def load_whitelist(file_path):
    """
    Load whitelist barcodes.

    Expected file format:
        - Each line contains a single barcode (e.g., "ATCGT").

    Parameters:
        file_path (str): Path to the whitelist file.

    Returns:
        list: A list of whitelist barcodes as strings.
    """
    with open(file_path, 'r') as f:
        return [line.strip() for line in f]


def load_observed_barcodes(file_path):
    """
    Load observed barcodes and their counts.

    Expected file format:
        - Tab-separated file.
        - Column 1: Observed barcode (e.g., "ATCGT").
        - Column 2: Count (integer).

    Example:
        ATCGT   100
        GCTAA   50

    Parameters:
        file_path (str): Path to the observed barcodes file.

    Returns:
        dict: A dictionary where keys are barcodes (str) and values are counts (int).
    """
    observed = {}
    with open(file_path, 'r') as f:
        for line in f:
            barcode, count = line.strip().split('\t')
            observed[barcode] = int(count)
    return observed


def load_quality_scores(file_path):
    """
    Load barcodes and their base quality scores.

    Expected file format:
        - Tab-separated file.
        - Column 1: Barcode (e.g., "ATCGT").
        - Columns 2+: Space-separated quality scores for each base (integers).

    Example:
        ATCGT   30 30 20 25 15
        GCTAA   25 25 30 30 20

    Parameters:
        file_path (str): Path to the quality scores file.

    Returns:
        dict: A dictionary where keys are barcodes (str) and values are lists of quality scores (list of int).
    """
    quality_scores = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            barcode = parts[0]
            scores = list(map(int, parts[1:]))
            quality_scores[barcode] = scores
    return quality_scores


def correct_barcodes(whitelist, observed, quality_scores, threshold=0.975):
    """
    Correct observed barcodes based on whitelist and posterior probabilities.

    Parameters:
        whitelist (list): List of whitelist barcodes.
        observed (dict): Dictionary of observed barcodes and their counts.
        quality_scores (dict): Dictionary of observed barcodes and their base quality scores.
        threshold (float): Posterior probability threshold for correction.

    Returns:
        dict: A dictionary of corrected barcodes and their counts.
    """
    corrected = {}

    # Count observed frequencies of whitelist barcodes
    whitelist_counts = {barcode: observed.get(barcode, 0) for barcode in whitelist}

    for obs_barcode, obs_count in observed.items():
        if obs_barcode in whitelist:
            # Direct match with whitelist
            corrected[obs_barcode] = corrected.get(obs_barcode, 0) + obs_count
            continue

        # Find whitelist barcodes within Hamming distance of 1
        candidates = [
            wl_barcode for wl_barcode in whitelist
            if hamming_distance(obs_barcode, wl_barcode) == 1
        ]

        posterior_probs = {}
        for wl_barcode in candidates:
            # Identify the differing position
            differing_positions = [
                i for i, (obs_char, wl_char) in enumerate(zip(obs_barcode, wl_barcode))
                if obs_char != wl_char
            ]
            if len(differing_positions) != 1:
                continue  # Skip invalid candidates

            differing_pos = differing_positions[0]
            # Get the quality score for the differing base
            quality_score = quality_scores.get(obs_barcode, [20] * len(obs_barcode))[differing_pos]
            error_prob = error_probability(quality_score)

            # Compute posterior probability
            prior_prob = whitelist_counts[wl_barcode]
            likelihood = error_prob
            posterior_probs[wl_barcode] = likelihood * prior_prob

        # Normalize posterior probabilities
        total_prob = sum(posterior_probs.values())
        for wl_barcode in posterior_probs:
            posterior_probs[wl_barcode] /= total_prob

        # Find the highest posterior probability
        if posterior_probs:
            best_candidate = max(posterior_probs, key=posterior_probs.get)
            best_prob = posterior_probs[best_candidate]

            if best_prob > threshold:
                # Replace observed barcode with whitelist barcode
                corrected[best_candidate] = corrected.get(best_candidate, 0) + obs_count
            else:
                # Keep the original barcode if no candidate passes the threshold
                corrected[obs_barcode] = corrected.get(obs_barcode, 0) + obs_count
        else:
            # No valid candidates, keep the original barcode
            corrected[obs_barcode] = corrected.get(obs_barcode, 0) + obs_count

    return corrected


def main():
    """
    Count the observed frequency of every barcode on the whitelist in the dataset.
    For every observed barcode in the dataset that is not on the whitelist and is at most one Hamming distance away from the whitelist sequences:
    Compute the posterior probability that the observed barcode did originate from the whitelist barcode but has a sequencing error at the differing base (by base quality score).
    Replace the observed barcode with the whitelist barcode that has the highest posterior probability (>0.975).

    Expected Input Files:
    1. Whitelist file (1st argument): Each line contains a barcode.
    2. Observed barcodes file (2nd argument): Tab-separated file with two columns:
        - Column 1: Observed barcode.
        - Column 2: Count of the barcode.
    3. Quality scores file (3rd argument): Tab-separated file with:
        - Column 1: Observed barcode.
        - Columns 2+: Quality scores for each base (integers).
    """
    parser = argparse.ArgumentParser(description="Correct barcodes using a whitelist and quality scores.")
    parser.add_argument("whitelist", help="Path to the whitelist file.")
    parser.add_argument("observed", help="Path to the observed barcodes file.")
    parser.add_argument("quality_scores", help="Path to the quality scores file.")
    parser.add_argument("output", help="Path to the output file.")
    parser.add_argument(
        "--threshold", type=float, default=0.975,
        help="Posterior probability threshold for correction (default: 0.975)."
    )
    args = parser.parse_args()

    # Load data
    whitelist = load_whitelist(args.whitelist)
    observed = load_observed_barcodes(args.observed)
    quality_scores = load_quality_scores(args.quality_scores)

    # Correct barcodes
    corrected_barcodes = correct_barcodes(whitelist, observed, quality_scores, args.threshold)

    # Save results
    with open(args.output, 'w') as f:
        for barcode, count in corrected_barcodes.items():
            f.write(f"{barcode}\t{count}\n")

    print(f"Corrected barcodes saved to {args.output}")


if __name__ == "__main__":
    main()

