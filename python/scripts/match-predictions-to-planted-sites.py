#
# Copyright John Reid 2011
#

"""
Matches predictions in MEME output format to details about planted sites in headers of fasta files.
"""

import logging
from stempy.planted_sites import parse_fasta_for_sites, parse_meme_output_for_sites, calculate_positives_and_negatives
from infpy.roc import RocCalculator

logging.basicConfig(level=logging.INFO)

if 3 != len(sys.argv):
    raise RuntimeError('USAGE: %s <fasta> <meme-output>' % sys.argv[0])
fasta_file = sys.argv[1]
meme_output = sys.argv[2]

records, true_sites = parse_fasta_for_sites(fasta_file)
predicted_sites = parse_meme_output_for_sites(meme_output)
TP, TN, FN, FP = calculate_positives_and_negatives(
    records, true_sites, predicted_sites)
roc = RocCalculator(TP, FP, TN, FN)
logging.info(roc)
