#
# Copyright John Reid 2012
#

"""
Test STEME on data from Hemant.
"""

from setup_environment import init_test_env, logging
init_test_env(__file__, level=logging.INFO)

import os, sys, stempy


#
# Run STEME
#
options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-run-steme-2')
options.bg_fasta_file = os.path.join(os.path.dirname(__file__), 'fasta', 'SOX2-t=20.fasta').encode(sys.stdin.encoding or 'ascii')
options.min_num_sites = options.max_num_sites = 10
options.bg_model_order = 4
options.use_seed = 'TCTGTGTTTGCATT'
options.prediction_Z_threshold = 1e-10
options.num_motifs = 1

algorithm = stempy.Algorithm(options)
fasta = os.path.join(os.path.dirname(__file__), 'fasta', 'PAX6-SOX2-t=20.fasta.masked')
algorithm(fasta.encode(sys.stdin.encoding or 'ascii'))
