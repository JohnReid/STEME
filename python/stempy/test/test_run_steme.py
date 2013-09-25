#
# Copyright John Reid 2011
#

"""
Test STEME on small-ish FASTA file.
"""

from setup_environment import init_test_env, logging
init_test_env(__file__, level=logging.INFO)

import os, sys, stempy


#
# Run STEME
#
options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-run-steme')
options.min_w = options.max_w = 8
options.write_em_stats = True
options.meme_like_output = 'meme.out'
algorithm = stempy.Algorithm(options)
fasta = os.path.join(os.path.dirname(__file__), 'fasta', 'T00759-tiny.fa')
algorithm(fasta.encode(sys.stdin.encoding or 'ascii'))
