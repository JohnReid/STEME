#
# Copyright John Reid 2012
#

"""
Test STEME when there are no site in the background.
"""

from setup_environment import init_test_env, logging
init_test_env(__file__, level=logging.INFO)

import os, sys, stempy


#
# Run STEME
#
options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-no-sites-in-bg')
options.min_w = 6
options.max_w = 14
options.bg_model_order = 3
options.bg_fasta_file = os.path.join(os.path.dirname(__file__), 'fasta', 'random-seqs-30-200.fasta').encode(sys.stdin.encoding or 'ascii')
fasta = os.path.join(os.path.dirname(__file__), 'fasta', 'dm01r.fasta')
algorithm = stempy.Algorithm(options)
algorithm(fasta.encode(sys.stdin.encoding or 'ascii'))
