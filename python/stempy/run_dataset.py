#
# Copyright John Reid 2009, 2010
#

"""
Code to run STEM/MEME algorithm on NTNU/Tompa data sets.
"""

import os
import logging
import stempy.tompa as tompa
import stempy.ntnu as ntnu
import stempy as stem
import stempy.meme as meme
from stempy import rationalise_predictions, ensure_dir_exists


def method_for_name(method_name):
    if 'STEM' == method_name:
        return stem
    elif 'MEME' == method_name:
        return meme
    else:
        raise RuntimeError('Unknown method(s): "%s"' % method_name)


def suite_for_name(suite_name):
    if 'NTNU' == suite_name:
        return ntnu
    elif 'Tompa' == suite_name:
        return tompa
    else:
        raise RuntimeError('Unknown benchmark suite: "%s"' % suite_name)


def output_dir_for_dataset(suite_name, data_set):
    "@return: A directory to output results of running dataset."
    return os.path.join('output', 'run_dataset', suite_name, *data_set)


def run_dataset(method_name, suite_name, data_set, fasta, options):
    "Run data set."
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    options.output_dir = output_dir_for_dataset(suite_name, data_set)
    ensure_dir_exists(options.output_dir)
    handler = logging.FileHandler(
        os.path.join(options.output_dir, '%s.log' % method_name))
    handler.setLevel(logging.INFO)
    handler.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)
    logging.info('%s is analysing data set: %s', method_name, data_set)
    predictions = method_for_name(method_name).Algorithm(options)(fasta)
    logging.info('%s predicted %d binding sites.',
                 method_name, len(predictions))
    predictions = list(rationalise_predictions(predictions))
    return predictions
