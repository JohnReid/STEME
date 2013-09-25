#
# Copyright John Reid 2009, 2010
#

"""
Code to run EM part of STEME algorithm.
"""

#
# Set up the logging
#
import logging
import sys
from cookbook.script_basics import setup_logging
setup_logging(__file__, level=logging.INFO)

# show_environment()

#
# Set up options
#
import stempy
options, args = stempy.parse_options(stempy.add_options)
if len(args) != 0:
    raise RuntimeError('USAGE: %s <options>', sys.argv[0])

W = 8
fasta_file = '/home/john/Data/GappedPssms/apr-2009/T99006trimRM.fa'
algorithm = stempy.Algorithm(options)
algorithm.initialise(fasta_file)
model = algorithm.create_default_model(W)
model.prior_num_sites = 18.276144706645898
model.lambda_ = 0.00037315491488373951
model.bs.pssm.log_probs.values()[:] = [
    [0.012967,  0.884511,  0.064057,  0.038465],
    [0.021795,  0.875048,  0.023177,  0.079979],
    [0.031394,  0.912065,  0.018154,  0.038387],
    [0.22118,  0.486244,  0.067231,  0.225346],
    [0.065432,  0.248463,  0.447243,  0.238861],
    [0.024532,  0.855213,  0.057096,  0.06316],
    [0.025442,  0.904468,  0.00916,  0.06093],
    [0.026273,  0.921693,  0.012526,  0.039508],
]
model.bs.recalculate()
model.bs.seed_pseudo_counts = options.em_seed_pseudo_counts / \
    4.  # divide by alphabet length just as MEME does.
EM = stempy.create_EM_descender(
    algorithm.data, model, options.epsilon, options.wnsites)
EM.using_sparse_Z = False
logging.info('Press <return>')
sys.stdin.readline()  # wait for user input
logging.info('EM.do_iteration()')
EM.do_iteration()
for window_start in xrange(algorithm.data.N - W + 1):
    z_sum = 0.
    for n in xrange(window_start, window_start + W):
        z = EM.get_Z(n)
        if z.first:
            z_sum += z.first.value
        if z.second:
            z_sum += z.second.value
    if z_sum > 1. + 1e-3:
        for n in xrange(window_start, window_start + W):
            z = EM.get_Z(n)
            if z.first:
                print '%d + %f' % (n, z.first.value)
            if z.second:
                print '%d - %f' % (n, z.second.value)
        raise ValueError(
            'Z have not been normalised. window start=%d; sum(z)=%f' %
            (window_start, z_sum))
