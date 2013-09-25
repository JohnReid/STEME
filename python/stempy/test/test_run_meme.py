#
# Copyright John Reid 2011
#

"""
Test STEME on small-ish FASTA file.
"""

from setup_environment import init_test_env, logging
init_test_env(__file__, level=logging.INFO)

import stempy, os
from stempy import meme

options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-run-meme')
fasta = os.path.join(os.path.dirname(__file__), 'fasta', 'T00759-tiny.fa')
meme_algorithm = meme.Algorithm(options)
meme_algorithm(fasta)
