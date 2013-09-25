#
# Copyright John Reid 2010, 2011
#

"""
Code to run the MEME EM algorithm.
"""

import logging
from cookbook.script_basics import setup_logging
setup_logging(__file__, level=logging.INFO)

from optparse import OptionParser
import stempy
import os
import sys
import tables
import math
from stempy import test_data


class MemeEM(tables.IsDescription):

    "A HDF5 class to represent the result of one EM run in MEME."
    dataset = tables.StringCol(32)
    cons0 = tables.StringCol(32)
    nsites0 = tables.Int32Col()
    niters = tables.Int32Col()
    em_time = tables.FloatCol()
    cons = tables.StringCol(32)
    nsites = tables.Int32Col()
    sig = tables.FloatCol()


def run_meme_on_start(fasta, seed, num_sites, score, options):
    import stempy.meme as meme
    _meme_cmd_args, stdoutdata, starts, _Zs, _thetas, _lambdas = meme.run_meme(fasta, options, extra_args=(
        '-nsites', str(num_sites), '-cons', seed, '-w', str(len(seed)))
    )
    if 1 != len(starts):
        raise RuntimeError(
            'Expecting only one start from MEME:\n%s' % stdoutdata)
    return starts[0]


def run(options):
    #
    # Create HDF5 table
    #
    h5file = tables.openFile(filename, mode="a", title="STEM/MEME data")
    try:
        meme_em_table = h5file.root.MEME.starts
    except tables.NoSuchNodeError:
        meme_group = h5file.createGroup(
            "/", 'MEME', 'Data about MEME runs')  # Create a new group under "/" (root)
        # Create one table on it
        meme_em_table = h5file.createTable(
            meme_group, 'starts', MemeEM, "Info on MEME starts.")

    #
    # Examine data sets
    #
    data_sets = test_data.data_sets_for_options(options)
    fasta_filenames = [test_data.fasta_filenames[ds] for ds in data_sets]
    for ds, fasta in zip(data_sets, fasta_filenames):
        num_seqs, num_bases = test_data.get_data_set_size(ds)
        logging.info(
            'Analysing data set: %16s; # seqs=%5d; # bases=%7d; %s', ds, num_seqs, num_bases, fasta)

    #
    # set up parallel stuff
    #
    from IPython.kernel import client
    tc = client.TaskClient()

    #
    # pass tasks to engines
    #
    logging.info('Passing tasks to engines')
    task_ids = []
    task_args = dict()
    task_data_set = dict()
    for data_set, fasta in zip(data_sets, fasta_filenames):
        for seed, num_sites, score in test_data.starts[data_set]:
            # only pass task if we don't have data in the table already
            if 0 == len(
                meme_em_table.getWhereList(
                    '(dataset=="%s") & (cons0=="%s") & (nsites0==%d)' % (
                        data_set, seed, math.trunc(num_sites)
                    )
                )
            ):
                options.output_dir = os.path.abspath(
                    os.path.join('output', 'meme-em', '%s-%03d' % (seed, num_sites)))
                os.path.exists(options.output_dir) or os.makedirs(
                    options.output_dir)
                args = (fasta, seed, num_sites, options)
                task = client.MapTask(
                    run_meme_on_start, (fasta, seed, num_sites, score, options))
                task_id = tc.run(task, block=False)
                task_ids.append(task_id)
                task_data_set[task_id] = data_set
                task_args[task_id] = args

    #
    # Get results from engines
    #
    logging.info('Blocking on %d results...', len(task_ids))
    meme_em = meme_em_table.row  # Fill the table with data
    for task_id in task_ids:
        start = tc.get_task_result(task_id, block=True)
        fasta, seed, num_sites, options = task_args[task_id]
        assert seed == start.cons0
        data_set = task_data_set[task_id]
        meme_em['dataset'] = data_set
        meme_em['cons0'] = start.cons0
        meme_em['nsites0'] = start.nsites0
        meme_em['niters'] = start.niters
        meme_em['em_time'] = start.em_time
        meme_em['cons'] = start.cons_after_em
        meme_em['nsites'] = start.nsites
        meme_em['sig'] = start.sig
        logging.info(
            '%s: cons0=%20s; nsites0=%3d; niters=%4d; elapsed=%7.1fs; per iteration=%6.2fs; cons=%20s; nsites0=%3d; sig=%e',
            data_set, start.cons0, start.nsites0, start.niters, start.em_time, start.em_time /
            start.niters, start.cons, start.nsites, start.sig
        )
        meme_em.append()
    h5file.close()  # Close (and flush) the HDF5 file


if '__main__' == __name__:
    #
    # parse options
    #
    parser = OptionParser()
    test_data.add_options(parser)
    stempy.add_options(parser)
    options, args = parser.parse_args()
    usage = 'USAGE: %s <options> <h5 file>' % sys.argv[0]
    if len(args) != 1:
        raise RuntimeError(usage)
    if len(args) < 1:
        raise RuntimeError(usage)
    filename = args[0]

    run(options)
