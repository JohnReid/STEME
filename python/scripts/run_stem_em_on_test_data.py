#
# Copyright John Reid 2010, 2011
#

"""
Code to test the STEME EM algorithm for stability with respect to the epsilon approximation parameter.
"""

import logging
from cookbook.script_basics import setup_logging
setup_logging(__file__, level=logging.INFO)

import stempy
import os
import sys
import tables
import math
from stempy import test_data


class Timings(tables.IsDescription):

    "A HDF5 class to represent timings."
    dataset = tables.StringCol(32)
    seed = tables.StringCol(32)
    nsites = tables.Int32Col()
    epsilon = tables.FloatCol()
    niters = tables.Int32Col()
    duration = tables.FloatCol()
    consensus = tables.StringCol(32)


class Dataset(tables.IsDescription):

    "A HDF5 class to represent data sets."
    dataset = tables.StringCol(32)
    size = tables.Int32Col()


def test_stability_and_speed(fasta, seed, num_sites, score, epsilon, options):
    import stempy
    import logging
    from cookbook.timer import Timer
    algorithm = stempy.Algorithm(options)
    algorithm.initialise(fasta)
    start = stempy.Start(
        seed=seed, num_sites=num_sites, score=score, model=None)
    options.epsilon = epsilon
    logging.debug(
        'Testing seed=%s; num_sites=%d epsilon=%f; data size=%d, fasta=%s',
        seed, num_sites, epsilon, algorithm.data.num_W_mers(len(seed)), fasta
    )
    algorithm.index.visit()  # pre-visit index to make sure is built correctly.
    with Timer(msg='run EM') as timer:
        em_result = algorithm.run_em_from_start(start)
    _duration = timer.duration
    _post_EM_consensus = stempy.consensus_from_pssm(
        em_result.model.bs.pssm.log_probs.values())
    return em_result.em_duration, em_result.cons_after_em, len(em_result.LLs)


hostname = os.uname()[1]
filename = "%s.h5" % hostname
epsilons = [0., .2, .4, .6, .8]
#epsilons = [0.,]


def run(options):
    #
    # Create HDF5 table
    #
    h5file = tables.openFile(filename, mode="a", title="Eta stability data")
    try:
        timings_table = h5file.root.EtaStability.timings
    except tables.NoSuchNodeError:
        group = h5file.createGroup(
            "/", 'EtaStability', 'Data about stability of epsilon')  # Create a new group under "/" (root)
        # Create one table on it
        timings_table = h5file.createTable(
            group, 'timings', Timings, "Timings table")

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
            options.max_num_sites = num_sites
            for epsilon in epsilons:
                # only pass task if we don't have data in the table already
                where_list = timings_table.getWhereList(
                    '(dataset=="%s") & (seed=="%s") & (nsites==%d) & (epsilon>%f) & (epsilon<%f)' % (
                        data_set, seed, math.trunc(
                            num_sites), epsilon - 1e-4, epsilon + 1e-4
                    )
                )
                if options.force or 0 == len(where_list):
                    # remove data if we already have it and are forcing new
                    # calculation
                    if len(where_list):
                        if 1 != len(where_list):
                            raise ValueError(
                                'Expecting to just find one existing row for this start.')
                        timings_table.removeRows(where_list[0])
                    options.output_dir = os.path.abspath(
                        os.path.join('output', 'epsilon-stability', '%s-%03d-%.1f' % (seed, num_sites, epsilon)))
                    os.path.exists(options.output_dir) or os.makedirs(
                        options.output_dir)
                    args = (fasta, seed, num_sites, epsilon, options)
                    # print fasta, seed, num_sites, score, epsilon, options
                    task = client.MapTask(test_stability_and_speed, (
                        fasta, seed, num_sites, score, epsilon, options))
                    task_id = tc.run(task, block=False)
                    task_ids.append(task_id)
                    task_data_set[task_id] = data_set
                    task_args[task_id] = args

    #
    # Get results from engines
    #
    logging.info('Blocking on %d results...', len(task_ids))
    timings = timings_table.row  # Fill the table with data
    for task_id in task_ids:
        duration, post_EM_consensus, num_iters = tc.get_task_result(
            task_id, block=True)
        fasta, seed, num_sites, epsilon, options = task_args[task_id]
        data_set = task_data_set[task_id]
        timings['dataset'] = data_set
        timings['seed'] = seed
        timings['nsites'] = num_sites
        timings['epsilon'] = epsilon
        timings['niters'] = num_iters
        timings['duration'] = duration
        timings['consensus'] = post_EM_consensus
        logging.info(
            '%20s; nsites=%3d; epsilon=%.1f; iters=%5d; elapsed=%7.1fs; per iteration=%6.2fs; %20s; %s',
            seed, num_sites, epsilon, num_iters, duration, duration /
            num_iters, post_EM_consensus, data_set
        )
        timings.append()
    h5file.close()  # Close (and flush) the HDF5 file


def time_per_iteration_per_base(timings):
    return timings['duration'] / timings['niters'] / test_data.get_num_w_mers(timings['dataset'], len(timings['seed']))


if '__main__' == __name__:

    def add_options(parser):
        parser.add_option("-f", "--force", action="store_true",
                          help="Always run even if already have results.")

    #
    # parse options
    #
    options, args = stempy.parse_options(
        lambda parser: (add_options(parser), stempy.add_options(
            parser), test_data.add_options(parser))
    )

    if len(args) != 0:
        raise RuntimeError('USAGE: %s <options>', sys.argv[0])

    run(options)

    #
    # Slice and dice data
    #
    h5file = tables.openFile(filename)
    timings_table = h5file.root.EtaStability.timings
    nsite_values = list(set(row['nsites'] for row in timings_table))
    nsite_values.sort()
    time_per_iteration_per_base = [
        [time_per_iteration_per_base(timings)
         for timings in timings_table.where('nsites == %d' % nsites)]
        for nsites in nsite_values
    ]
    h5file.close()  # Close the HDF5 file
