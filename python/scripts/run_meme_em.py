#
# Copyright John Reid 2010
#

"""
Code to run the MEME EM algorithm.
"""

import logging
from cookbook.script_basics import setup_logging
setup_logging(__file__, level=logging.INFO)

import stempy
import stempy.meme as meme
import sys


#
# parse options
#
options, args = stempy.parse_options(stempy.add_options)
if len(args) != 3:
    raise RuntimeError('USAGE: %s <options> fasta seed num_sites', sys.argv[0])
fasta = args.pop(0)
seed = args.pop(0)
num_sites = int(args.pop(0))

meme_cmd_args, stdoutdata, starts, Zs = meme.run_meme(
    fasta,
    options,
    extra_args=(
        '-nsites', str(num_sites), '-cons', seed, '-w', str(
            len(seed)), "-print_z"
    )
)
