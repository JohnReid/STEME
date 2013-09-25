#
# Copyright John Reid 2009, 2010, 2011
#

"""
Check STEME EM against MEME EM.
"""

#
# Set up the logging
#
import logging
from cookbook.script_basics import setup_logging
setup_logging(__file__, level=logging.INFO)

import os
import tables
import numpy as np
import pylab as pl
import sys
from stempy import test_data
from optparse import OptionParser
from collections import defaultdict
from cookbook.pylab_utils import violin_plot

from setup_pylab import setup_pylab_for_tex, correct_axes
setup_pylab_for_tex()


def get_meme_row_for_stem_row(stem_row):
    "@return: The MEME row corresponding to the STEME row."
    condition = '(dataset == "%s") & (cons0 == "%s") & (nsites0 == %d)' % (
        stem_row['dataset'], stem_row['seed'], stem_row['nsites']
    )
    meme_rows = meme_em_table.getWhereList(condition)
    if not meme_rows:
        return None
    elif len(meme_rows) > 1:
        raise RuntimeError('More than one MEME row matches the STEME row.')
    else:
        return meme_em_table[meme_rows[0]]


def get_epsilon_0_row_for_stem_row(stem_row):
    "@return: The STEME row with epsilon=0 corresponding to the STEME row."
    condition = '(dataset == "%s") & (seed == "%s") & (nsites == %d) & (epsilon < 0.01)' % (
        stem_row['dataset'], stem_row['seed'], stem_row['nsites']
    )
    rows = stem_em_table.getWhereList(condition)
    if len(rows) == 0:
        return None
    elif len(rows) > 1:
        raise RuntimeError(
            'More than one STEME row with epsilon=0 matches the STEME row.')
    else:
        return stem_em_table[rows[0]]


def epsilon_to_index(epsilon):
    return int(epsilon * 1000)


def index_to_epsilon(index):
    return float(index) / 1000.


def hamming_distance(s1, s2):
    "@return: The Hamming distance between s1 and s2."
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def output_filename(name):
    return os.path.join(output_dir, name)


parser = OptionParser()
test_data.add_options(parser)
options, args = parser.parse_args()
if len(args) < 1:
    raise RuntimeError('USAGE: %s <h5 file>' % sys.argv[0])
data_sets = test_data.data_sets_for_options(options)
h5_filename = args[0]
if not data_sets:
    raise ValueError('No data sets specified in options.')

output_dir = os.path.join('output', 'STEM-vs-MEME')
os.path.exists(output_dir) or os.makedirs(output_dir)

h5file = tables.openFile(h5_filename)

#
# Get HDF5 tables
#
meme_em_table = h5file.root.MEME.starts
stem_em_table = h5file.root.EtaStability.timings

default_epsilon = .4
default_epsilon_index = epsilon_to_index(default_epsilon)

epsilon_mismatches = defaultdict(list)          # Seeds that have mismatches
epsilon_hamming = defaultdict(list)             # Hamming distances
# Fraction of consensus that has a mismatch
epsilon_fraction_mismatch = defaultdict(list)
# Seeds that have mismatches to epsilon=0
epsilon_0_mismatches = defaultdict(list)
# Fraction of consensus that has a mismatch to epsilon=0
epsilon_0_fraction_mismatch = defaultdict(list)
epsilon_iter_rel_speed = defaultdict(list)
epsilon_rel_speed = defaultdict(list)
width_iter_rel_speed = defaultdict(list)
num_sites_iter_rel_speed = defaultdict(list)
stem_runtime_by_size = defaultdict(list)
meme_runtime_by_size = defaultdict(list)
stem_itertime_by_size = defaultdict(list)
meme_itertime_by_size = defaultdict(list)
stem_num_iters = list()
meme_num_iters = list()

for stem_row in stem_em_table:
    epsilon_index = epsilon_to_index(stem_row['epsilon'])
    data_set = stem_row['dataset']
    if data_set in data_sets:
        stem_consensus = stem_row['consensus']
        W = len(stem_consensus)
        meme_row = get_meme_row_for_stem_row(stem_row)
        # print 'STEM:', stem_row
        # print 'MEME:', meme_row

        epsilon_0_row = get_epsilon_0_row_for_stem_row(stem_row)
        hamming_0 = hamming_distance(
            stem_consensus, epsilon_0_row['consensus'])
        epsilon_0_mismatches[epsilon_index].append(hamming_0 > 0)
        epsilon_0_fraction_mismatch[epsilon_index].append(hamming_0 / float(W))

        if meme_row:
            meme_consensus = meme_row['cons']
            hamming = hamming_distance(meme_row['cons'], stem_row['consensus'])
            # print meme_consensus, stem_consensus, stem_row['dataset'], hamming
            # print stem_row
# if hamming > 0 and epsilon_index == 0: # print details of mismatches when epsilon=0
#                print stem_row
#                print meme_consensus
#                print stem_consensus
            epsilon_hamming[epsilon_index].append(hamming)
            epsilon_fraction_mismatch[epsilon_index].append(
                hamming / float(len(meme_row['cons'])))
            epsilon_mismatches[epsilon_index].append(hamming > 0)
            rel_speed = stem_row['duration'] / meme_row['em_time']
            iter_rel_speed = rel_speed / \
                stem_row['niters'] * meme_row['niters']
            stem_num_iters.append(stem_row['niters'])
            meme_num_iters.append(meme_row['niters'])
            epsilon_rel_speed[epsilon_index].append(rel_speed)
            epsilon_iter_rel_speed[epsilon_index].append(iter_rel_speed)
            if epsilon_index == default_epsilon_index:
                width_iter_rel_speed[W].append(iter_rel_speed)
                num_sites_iter_rel_speed[
                    meme_row['nsites']].append(iter_rel_speed)
                num_seqs, num_bases = test_data.get_data_set_size(data_set)
                stem_runtime_by_size[num_bases].append(
                    np.log10(stem_row['duration']))
                meme_runtime_by_size[num_bases].append(
                    np.log10(meme_row['em_time']))
                stem_itertime_by_size[num_bases].append(
                    np.log10(stem_row['duration']) - np.log10(stem_row['niters']))
                meme_itertime_by_size[num_bases].append(
                    np.log10(meme_row['em_time']) - np.log10(meme_row['niters']))


epsilon_indices = epsilon_mismatches.keys()
epsilon_indices.sort()
epsilon_range = np.arange(len(epsilon_indices))
epsilons = map(index_to_epsilon, epsilon_indices)
str_epsilons = ['%.1f' % e for e in epsilons]

Ws = width_iter_rel_speed.keys()
Ws.sort()
W_range = np.arange(len(Ws))
str_Ws = map(str, Ws)
num_sites_values = num_sites_iter_rel_speed.keys()
num_sites_values.sort()
num_sites_range = np.arange(len(num_sites_values))


boxplots = False
pl.rcParams['patch.facecolor'] = 'grey'
pl.close('all')


#
# Plot the runtimes by size
#
pl.figure()
correct_axes()
dataset_sizes = stem_runtime_by_size.keys()
dataset_sizes.sort()
meme_values = [meme_runtime_by_size[k] for k in dataset_sizes]
meme_means = map(np.mean, meme_values)
meme_stddevs = map(np.std, meme_values)
pl.errorbar(np.log10(dataset_sizes), meme_means,
            meme_stddevs, label='MEME', fmt='r:')
stem_values = [stem_runtime_by_size[k] for k in dataset_sizes]
stem_means = map(np.mean, stem_values)
stem_stddevs = map(np.std, stem_values)
pl.errorbar(np.log10(dataset_sizes), stem_means,
            stem_stddevs, label='STEME', color='k')
pl.ylabel('\(\log_{10} secs/iteration\)')
pl.xlabel('\(\\log_{10} \# bases\)')
pl.legend(loc='lower right')
pl.savefig(output_filename('runtimes-by-size.eps'))
pl.close()


#
# Plot the itertimes by size
#
pl.figure()
correct_axes()
dataset_sizes = stem_itertime_by_size.keys()
dataset_sizes.sort()
meme_values = [meme_itertime_by_size[k] for k in dataset_sizes]
meme_means = map(np.mean, meme_values)
meme_stddevs = map(np.std, meme_values)
pl.errorbar(np.log10(dataset_sizes), meme_means,
            meme_stddevs, label='MEME', fmt='r:')
stem_values = [stem_itertime_by_size[k] for k in dataset_sizes]
stem_means = map(np.mean, stem_values)
stem_stddevs = map(np.std, stem_values)
pl.errorbar(np.log10(dataset_sizes), stem_means,
            stem_stddevs, label='STEME', color='black')
pl.ylabel('\(\log_{10}\) seconds')
pl.xlabel('\(\\log_{10}\) \# bases')
pl.legend(loc='lower right')
pl.savefig(output_filename('itertimes-by-size.eps'))
pl.close()


#
# Plot the epsilon mismatch rates
#
pl.figure()
correct_axes()
values = [sum(epsilon_mismatches[epsilon_index]) / float(len(epsilon_mismatches[epsilon_index]))
          for epsilon_index in epsilon_indices]
pl.bar(epsilon_range - .4, values, alpha=.5)
#pl.ylim(0., 1.)
pl.xlim(-.5, -.5 + len(epsilon_range))
pl.ylabel('mismatch rate')
pl.xlabel('\(\epsilon\)')
ax = pl.gca()
ax.set_xticks(epsilon_range)
ax.set_xticklabels(str_epsilons)
pl.savefig(output_filename('epsilon-mismatch.eps'))
pl.close()

#
# Plot the epsilon fraction mismatch
#
pl.figure()
correct_axes()
values = [np.mean(epsilon_fraction_mismatch[epsilon_index])
          for epsilon_index in epsilon_indices]
#std = [np.std(epsilon_fraction_mismatch[epsilon_index]) for epsilon_index in epsilon_indices]
pl.bar(epsilon_range - .4, values, ecolor='grey', alpha=.5)
#pl.ylim(0., 1.)
pl.xlim(-.5, -.5 + len(epsilon_range))
pl.ylabel('mismatch fraction')
pl.xlabel('\(\epsilon\)')
ax = pl.gca()
ax.set_xticks(epsilon_range)
ax.set_xticklabels(str_epsilons)
pl.savefig(output_filename('epsilon-fraction-mismatch.eps'))
pl.close()

#
# Plot the epsilon 0 mismatch rates
#
pl.figure()
correct_axes()
values = [sum(epsilon_0_mismatches[epsilon_index]) / float(len(epsilon_mismatches[epsilon_index]))
          for epsilon_index in epsilon_indices]
pl.bar(epsilon_range - .4, values, alpha=.5)
#pl.ylim(0., 1.)
pl.xlim(-.5, -.5 + len(epsilon_range))
pl.ylabel('mismatch rate')
pl.xlabel('\(\epsilon\)')
ax = pl.gca()
ax.set_xticks(epsilon_range)
ax.set_xticklabels(str_epsilons)
pl.savefig(output_filename('epsilon-0-mismatch.eps'))
pl.close()

#
# Plot the epsilon 0 fraction mismatch
#
pl.figure()
correct_axes()
values = [np.mean(epsilon_0_fraction_mismatch[epsilon_index])
          for epsilon_index in epsilon_indices]
#std = [np.std(epsilon_fraction_mismatch[epsilon_index]) for epsilon_index in epsilon_indices]
pl.bar(epsilon_range - .4, values, ecolor='grey', alpha=.5)
pl.xlim(-.5, -.5 + len(epsilon_range))
pl.ylabel('mismatch fraction')
pl.xlabel('\(\epsilon\)')
ax = pl.gca()
ax.set_xticks(epsilon_range)
ax.set_xticklabels(str_epsilons)
pl.savefig(output_filename('epsilon-0-fraction-mismatch.eps'))
pl.close()

#
# Plot the epsilon iteration relative speeds in a violin plot
#
pl.figure()
correct_axes()
violin_plot(
    pl.gca(),
    [np.log10(epsilon_iter_rel_speed[epsilon_index])
     for epsilon_index in epsilon_indices],
    epsilon_range,
    bp=boxplots,
    facecolor='grey',
    alpha=1
)
pl.ylabel('\(\log_{10}\) relative speed')
pl.xlabel('\(\epsilon\)')
ax = pl.gca()
ax.set_xticks(epsilon_range)
ax.set_xticklabels(str_epsilons)
pl.savefig(output_filename('epsilon-iter-rel-speed.eps'))
pl.close()

#
# Plot the width iteration relative speeds in a violin plot
#
pl.figure()
correct_axes()
violin_plot(
    pl.gca(),
    [np.log10(width_iter_rel_speed[key]) for key in Ws],
    W_range,
    bp=boxplots,
    facecolor='grey',
    alpha=1
)
pl.ylabel('\(\log_{10}\) relative speed')
pl.xlabel('W')
ax = pl.gca()
ax.set_xticks(W_range)
ax.set_xticklabels(str_Ws)
pl.savefig(output_filename('width-iter-rel-speed.eps'))
pl.close()

#
# Plot the num-sites iteration relative speeds in a violin plot
#
pl.figure()
correct_axes()
violin_plot(
    pl.gca(),
    [np.log10(num_sites_iter_rel_speed[key]) for key in num_sites_values],
    num_sites_range,
    bp=boxplots,
    facecolor='grey',
    alpha=1
)
pl.ylabel('iteration log relative speed')
pl.xlabel('\\# sites')
ax = pl.gca()
ax.set_xticks(num_sites_range)
ax.set_xticklabels(map(str, num_sites_values))
pl.savefig(output_filename('num-sites-iter-rel-speed.eps'))
pl.close()

#
# Plot the epsilon relative speeds in a violin plot
#
pl.figure()
correct_axes()
violin_plot(
    pl.gca(),
    [np.log10(epsilon_rel_speed[epsilon_index])
     for epsilon_index in epsilon_indices],
    epsilon_range,
    bp=boxplots,
    facecolor='grey',
    alpha=1
)
pl.ylabel('log relative speed')
pl.xlabel('\(\epsilon\)')
ax = pl.gca()
ax.set_xticks(epsilon_range)
ax.set_xticklabels(str_epsilons)
pl.savefig(output_filename('epsilon-rel-speed.eps'))
pl.close()

# pl.show()
