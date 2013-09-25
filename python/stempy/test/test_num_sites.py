#
# Copyright John Reid 2011
#

"""
Test STEME gets the number of sites correct.
"""

from setup_environment import init_test_env, logging, fasta_dir
init_test_env(__file__, level=logging.INFO)


import stempy, os
from stempy.planted_sites import parse_meme_output_for_sites

#
# Set up the options
#
site = 'AAGGTTCCTTGGAATT'
W = len(site)
fasta_file = os.path.join(fasta_dir(), 'random-seqs-4-sites.fasta')
options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-num-sites')
options.bg_model_order = 0
options.min_w = options.max_w = W
options.min_num_sites = 2
options.max_num_sites = 10
options.meme_like_output = 'test-num-sites.txt'
meme_output = os.path.join(options.output_dir, options.meme_like_output)


#
# Run the STEME algorithm
#
algorithm = stempy.Algorithm(options)
algorithm(fasta_file)

#
# Make sure we choose a motif that predicted 4 sites
#    
predicted_sites = parse_meme_output_for_sites(meme_output)
for seq, sites in predicted_sites.iteritems():
    for i, _id in enumerate(algorithm.input_sequences.ids):
        if _id.startswith(seq):
            break
    for site in sites:
        global_pos = algorithm.input_sequences.data.pos_globalise(i, site.first)
        logging.info('%2d %3d %s %s', i, site.first, algorithm.input_sequences.data.get_W_mer(W, global_pos), seq)
assert 4 == len(predicted_sites)

