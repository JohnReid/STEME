#
# Copyright John Reid 2009, 2010
#

"""
Code to check the Z's that STEME and MEME calculate in the EM algorithm.
"""

#
# Set up the logging
#
import logging
import os
import sys
from cookbook.script_basics import setup_logging
setup_logging(__file__, level=logging.INFO)
import stempy
import stempy.meme as meme
import numpy as np

# show_environment()

#
# Set up options
#
options, args = stempy.parse_options(stempy.add_options)
options.output_dir = os.path.abspath(
    os.path.join('output', 'STEM-vs-MEME', 'Z'))
if len(args) != 3:
    raise RuntimeError('USAGE: %s <options> fasta seed num_sites', sys.argv[0])
fasta_file = args.pop(0)
seed = args.pop(0)
num_sites = int(args.pop(0))
if options.epsilon > 0.:
    raise ValueError('No point running this script with epsilon > 0.')


#
# Run MEME
#
meme_cmd_args, meme_stdoutdata, meme_starts, meme_Zs, meme_thetas, meme_lambdas = meme.run_meme(
    fasta_file,
    options,
    extra_args=(
        '-nsites', str(num_sites), '-cons', seed, '-w', str(
            len(seed)), '-print_z', '-print_all'
    )
)
options.max_iters = len(meme_Zs)


def close(z1, z2, rtol=1e-2, atol=1e-1):
    return abs(z1 - z2) <= (atol + rtol * abs(z2))


def stem_cb(**kwargs):
    "Callback for STEME EM."
    stem = kwargs['stem']
    EM = kwargs['em']
    _iter = kwargs['iter']
    stem_theta = kwargs['theta']
    if iter < len(meme_Zs):
        meme_Z = meme_Zs[iter]
        meme_theta = np.asarray(meme_thetas[iter][len(seed):])
        # print np.log(meme_theta[-1,-1])
        # print meme_theta
        meme_lambda = meme_lambdas[iter]
        euclidean = np.sqrt(
            ((stem_theta - meme_theta) ** 2).sum()) / stem_theta.size
        lambda_diff = EM.model.lambda_ - meme_lambda
        logging.info(
            'STEME iteration %d: theta distance=%f; lambda abs. diff.=%f; ; lambda rel. diff.=%f',
            iter, euclidean, lambda_diff, lambda_diff / meme_lambda)
        for (seq, pos, strand), meme_z in meme_Z.iteritems():
            logging.debug('%5d seq; %7d pos; %s', seq, pos, strand)
            global_pos = stem.data.pos_globalise(seq, pos)
            stem_z_pair = EM.get_Z(global_pos)
            if '+' == strand:
                stem_z = stem_z_pair.first
            else:
                stem_z = stem_z_pair.second
            if stem_z:
                stem_z = stem_z.value
                if not close(meme_z, stem_z):
                    W_mer = stem.data.get_W_mer(8, global_pos)
                    raise ValueError(
                        'Zs do not match at (%d, %d, %s) [%d]; %s; MEME %f != %f STEM; abs. diff.=%e; rel. diff.=%e' % (
                            seq, pos, strand, global_pos, W_mer, meme_z, stem_z, meme_z -
                            stem_z, (meme_z - stem_z) / stem_z
                        )
                    )
    else:
        logging.info('No MEME results for iteration %d', iter)

#
# Run STEM
#
algorithm = stempy.Algorithm(options)
algorithm.initialise(fasta_file)
start = stempy.Start(num_sites=num_sites, score=0., seed=seed, model=None)
EM, _expected_sites, LLs, model = algorithm.run_em_from_start(
    start, callback=stem_cb)
