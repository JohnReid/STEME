#
# Copyright John Reid 2009, 2010
#

"""
Code to run STEM/MEME algorithm on NTNU/Tompa data sets.
"""

#
# Set up the logging
#
from cookbook.script_basics import setup_logging
import logging
setup_logging(__file__, level=logging.INFO)

import sys
import stempy as stem
import stempy.meme as meme
from stempy.run_dataset import method_for_name, suite_for_name


if '__main__' == __name__:
    #
    # get the method from the command line
    #
    if len(sys.argv) < 2:
        raise RuntimeError('No method(s) specified')
    method_names = sys.argv[1].split(':')
    sys.argv.pop(1)
    map(method_for_name, method_names)

    #
    # get the suite from the command line
    #
    if len(sys.argv) < 2:
        raise RuntimeError('No benchmark suite specified')
    suite_names = sys.argv[1].split(':')
    sys.argv.pop(1)
    map(suite_for_name, suite_names)

    #
    # Parse the options
    #
    from optparse import OptionParser
    option_parser = OptionParser()
    option_parser.add_option(
        "--num-threads",
        dest="num_threads",
        default=3,
        type='int',
        help="Number of threads to run jobs on."
    )
    option_parser.add_option("--data-sets", action="append")
    stem.add_options(option_parser)
    meme.add_options(option_parser)
    options = stem.parse_options(option_parser=option_parser)
    stem.turn_on_google_profiling_if_asked_for(options)

    # for each method and suite
    for method_name in method_names:
        for suite_name in suite_names:

            suite = suite_for_name(suite_name)
            method = method_for_name(method_name)

            predictions_by_dataset = []
            import cookbook.function_as_task as F

            def do_work(task):
                method_name, suite_name, data_set, options = task
                logging.info('Running %s', task)
                suite = suite_for_name(suite_name)
                fasta = suite.get_fasta_file_for_data_set(data_set)
                predictions = F.run_as_subprocess(
                    "stempy.run_dataset",
                    "run_dataset",
                    method_name,
                    suite_name,
                    data_set,
                    fasta,
                    options
                )
                predictions_by_dataset.append((data_set, fasta, predictions))

            q = F.create_queue(options.num_threads, do_work)

            #
            # Run the method on the data set(s)
            #
            method = method.Algorithm(options)

            # if any data sets specified on command line, just run those
            if options.data_sets:
                for ds in options.data_sets:
                    data_set = ds.split(':')
                    task = method_name, suite_name, data_set, options
                    q.put(task)

            # run all data sets
            else:
                for t in suite.get_data_set_types():
                    for data_set in suite.get_data_sets_for_type(t):
                        task = method_name, suite_name, data_set, options
                        q.put(task)

            # block until all tasks are done
            q.join()

            # write the predictions out...
            results_formatter = suite.ResultsFormatter(method_name)
            for data_set, fasta, predictions in predictions_by_dataset:
                results_formatter.start_data_set(data_set, fasta)
                for seq, interval in predictions:
                    if interval.empty():
                        logging.warning(
                            '%s returned empty interval: %s %s', method_name, seq, interval)
                    else:
                        results_formatter.write_prediction(seq, interval)
            results_formatter.finalise()

    stem.turn_off_google_profiling_if_asked_for(options)
