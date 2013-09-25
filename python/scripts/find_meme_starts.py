#
# Copyright John Reid 2010
#

"""
Code to find MEME starts for various data sets.
"""

import logging
from cookbook.script_basics import setup_logging
setup_logging(__file__, level=logging.INFO)

import stempy.meme
import os
import sys
from stempy import test_data

# parse options
options, args = stempy.parse_options(stempy.add_options)
if len(args) != 1:
    raise RuntimeError(
        'USAGE: %s <options> data_sets_size=[tiny|small|large]', sys.argv[0])
data_sets_size = args[0]

# work out which data sets we will run on
data_sets = getattr(test_data, '%s_data_sets' % data_sets_size)
fasta_filenames = [test_data.fasta_filenames[ds] for ds in data_sets]
for ds, fasta in zip(data_sets, fasta_filenames):
    logging.info('Analysing data set: %16s: %s', ds, fasta)


def get_starts(data_set, fasta):
    "Get the starts for the data set."
    stat_info = os.stat(fasta)
    logging.info('Finding starts in %s; size=%d', fasta, stat_info.st_size)
    options.output_dir = os.path.abspath(
        os.path.join('output', 'STEM-MEME', data_set))
    os.path.exists(options.output_dir) or os.makedirs(options.output_dir)
    _eme_cmd_args, _stdoutdata, starts, _Zs, _thetas, _lambdas = stempy.meme.run_meme(
        fasta, options)
    return starts

# do the parallel stuff
from IPython.kernel import client
mec = client.MultiEngineClient()
mec.push(dict(options=options))
tc = client.TaskClient()
results = tc.map(get_starts, data_sets, fasta_filenames)

# print results
for ds, starts in zip(data_sets, results):
    print "    '%s' : [" % ds
    print '\n'.join(
        "        ('%s', %.0f, %.2f)," % (start.cons0, start.nsites0, start.sig)
        for start in starts
    )
    print "    ],"
