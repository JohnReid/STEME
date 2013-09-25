#
# Copyright John Reid 2011
#

"""
Test STEME's performance on some small data sets of random sequences with planted sites.
"""

from setup_environment import init_test_env, logging, fasta_dir
init_test_env(__file__, level=logging.INFO)

import stempy, os
from stempy.planted_sites import parse_meme_output_for_sites
#from infpy.roc import RocCalculator
#from optparse import OptionParser

#rocs = dict()
#meme_rocs = dict()

#
# Set up the options
#
fasta_file = os.path.join(fasta_dir(), 'random-seqs-two-motifs.fasta')
options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-2-motifs')
options.min_w = 8
options.max_w = 10
options.num_motifs = 2
options.meme_like_output = 'two-motif-test-meme.txt'
meme_output = os.path.join(options.output_dir, options.meme_like_output)


#
# Run the STEME algorithm
#
algorithm = stempy.Algorithm(options)
algorithm(fasta_file)

#
# Make sure we can parse output with 2 motifs in it
#    
predicted_sites = parse_meme_output_for_sites(meme_output)

#
# Calculate the consensuses for the 2 motifs
#
consensuses = [
    stempy.consensus_from_pssm(motif.model.bs.pssm.log_probs.values())
    for motif in algorithm.motifs
]
assert consensuses[0] == 'AAACTCACTC' or stempy.reverse_complement(consensuses[0]) == 'AAACTCACTC', consensuses[0]
assert consensuses[1] == 'AACCTGTG'   or stempy.reverse_complement(consensuses[1]) == 'AACCTGTG', consensuses[1]
