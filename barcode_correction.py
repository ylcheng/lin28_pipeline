"""
Count the observed frequency of every barcode on the whitelist in the dataset.
For every observed barcode in the dataset that is not on the whitelist and is at most one Hamming distance away from the whitelist sequences:
Compute the posterior probability that the observed barcode did originate from the whitelist barcode but has a sequencing error at the differing base (by base quality score).
Replace the observed barcode with the whitelist barcode that has the highest posterior probability (>0.975).
"""
