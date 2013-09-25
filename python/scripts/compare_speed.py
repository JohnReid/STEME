#
# Copyright John Reid 2010
#

"""
Compare the speed of MEME and STEM.
"""


#
# Set up the logging and options
#
import logging
import os
import time
import pylab as P
from cookbook.dicts import DictOf
import cookbook.pylab_utils as pylab_utils
from cookbook.script_basics import setup_logging
setup_logging(__file__, level=logging.INFO)

import stempy
import stempy as stem
import stempy.meme as meme
stempy.Pssm.max_samples_in_E_value_calc = 1000


def add_options(option_parser):
    stem.add_options(option_parser)
    meme.add_options(option_parser)
options, args = stempy.parse_options(add_options)


def small_fastas():
    "@return: A few small fasta filenames."
    return [
        os.path.abspath(
            os.path.join(os.path.dirname(__file__), '../../fasta/%s.fa')) % f
        for f in (
            'T00759-tiny',
            'T00759-small',
            'T00759trimRM-test-x2',
        )
    ]


def big_fastas():
    "@return: A few big fasta files."
    return [
        '/home/john/Data/GappedPssms/apr-2009/T99006trimUN.fa',
        '/home/john/Data/GappedPssms/apr-2009/T99005trimUN.fa',
        '/home/john/Data/GappedPssms/apr-2009/T00759trimUN.fa',
        '/home/john/Data/GappedPssms/apr-2009/T00671trimUN.fa',
        '/home/john/Data/GappedPssms/Full-Sp1/Sp1-1000000.fa',
    ]


# fasta = '/home/john/Data/NTNU-TF-search-dataset/datasets/model_real/M00724.fas'
# fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../fasta/T00759trimRM-test-x2.fa'))
# fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../fasta/T00759-tiny.fa'))
output_dir = os.path.abspath(os.path.join('output', 'speed-test'))


def timeit(fn):
    "Time the execution of a function."
    start_time = time.time()
    fn()
    return time.time() - start_time


def time_method(method):
    "Time a method/width pair."
    logging.info('Timing method %s for width %d', method, W)
    return timeit(lambda: method.Algorithm(options)(fasta))

Ws = [6, 8, 10, 12]
#Ws = [6, 8]
#fastas = small_fastas()
fastas = big_fastas()


def save_timings():
    "Graph the timings."
    P.figure()
    for W, colour in zip(Ws, pylab_utils.simple_colours):
        P.loglog(fasta_sizes, stem_timings[
                 W], label='STEME W=%d' % W, ls='-', color=colour)
        P.loglog(fasta_sizes, meme_timings[
                 W], label='MEME W=%d' % W, ls='-.', color=colour)
        P.legend(loc='upper left')
        P.xlabel('\\# bases in data set')
        P.ylabel('seconds')
    P.savefig(os.path.join(output_dir, 'timings.eps'))
    P.savefig(os.path.join(output_dir, 'timings.png'))
    P.close()


#
# do the timings.
#
pylab_utils.set_rcParams_for_latex()
stem_timings = DictOf(list)
meme_timings = DictOf(list)
fasta_sizes = []
for fasta in fastas:
    for W in Ws:
        options.min_w = options.max_w = W
        options.output_dir = os.path.join(
            output_dir, 'W=%02d-%s' % (W, stempy.basename_wo_ext(fasta)))
        stempy.ensure_dir_exists(options.output_dir)
        stem_algorithm = stem.Algorithm(options)
        meme_algorithm = meme.Algorithm(options)
        stem_timings[W].append(timeit(lambda: stem_algorithm(fasta)))
        meme_timings[W].append(timeit(lambda: meme_algorithm(fasta)))
    fasta_sizes.append(stem_algorithm.num_bases)
    save_timings()
