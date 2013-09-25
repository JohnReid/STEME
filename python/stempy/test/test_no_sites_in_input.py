#
# Copyright John Reid 2012
#

"""
Test STEME when there are no sites in the input.
"""

from setup_environment import init_test_env, logging
init_test_env(__file__, level=logging.INFO)

import os, sys, stempy


#
# Run STEME
#
options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-no-sites-in-input')
options.min_w = 6
options.max_w = 6
options.num_motifs = 5
options.bg_model_order = 3
fasta = os.path.join(os.path.dirname(__file__), 'fasta', 'no-sites-in-input.fasta')
algorithm = stempy.Algorithm(options)
algorithm(fasta.encode(sys.stdin.encoding or 'ascii'))
