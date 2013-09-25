#
# Copyright John Reid 2009, 2010
#

"""
Package to handle NTNU benchmark data/test sets.

U{http://tare.medisin.ntnu.no/pages/tools.php}
"""

import os
import glob
import sys
import tempfile
import tarfile
import shutil
from itertools import imap
from stempy import _read_fasta

_data_dir = '/home/john/Data/NTNU-TF-search-dataset'


def get_data_set_types():
    "@return: The types of data sets."
    return ('algorithm_markov', 'algorithm_real', 'model_real')


def data_sets_dir():
    "@return: The data sets directory."
    return os.path.join(_data_dir, 'datasets')


def data_sets_dir_for_type(type):
    "@return: The data sets directory for the type."
    return os.path.join(data_sets_dir(), type)


def test_sets_dir():
    "@return: The test sets directory."
    return os.path.join(_data_dir, 'testsets')


def test_sets_dir_for_type(type):
    "@return: The test sets directory for the type."
    return os.path.join(test_sets_dir(), type)


def data_set_from_fasta(fasta):
    "@return: The data set name given the fasta filename."
    dir, file = os.path.split(fasta)
    dir, type = os.path.split(dir)
    return type, os.path.splitext(file)[0]


def get_data_sets_for_type(type):
    "@return: A list of the data sets for the given type."
    return map(data_set_from_fasta, glob.glob(os.path.join(data_sets_dir_for_type(type), '*.fas')))


def get_fasta_file_for_data_set(dataset):
    "@return: The fasta filename for the given data set."
    type, name = dataset
    return os.path.join(data_sets_dir_for_type(type), '%s.fas' % name)


def get_test_file_for_data_set(type, dataset):
    "@return: The predictions filename for the given data set."
    return os.path.join(test_sets_dir_for_type(type), '%s_pred.txt' % dataset)


def ensure_dir_exists(directory):
    "Make sure a directory exists."
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        if not os.path.isdir(directory):
            raise RuntimeError(
                '%s already exists and is not a directory.' % directory)


class ResultsFormatter(object):

    "Formats binding site predictions into a format suitable for submitting to the web service."

    def __init__(self, method_name):
        "Initialise."
        self.file = 'NTNU-%s-results.tar.gz' % method_name
        self.tmpdir = tempfile.mkdtemp()
        self.result_files = []
        self.f = None

    def start_data_set(self, dataset, fasta):
        "Called when starting to format results of particular dataset."
        self.seq_dict, self.seq_indices, self.seq_names = _read_fasta(fasta)
        type, name = dataset
        dir = os.path.join(self.tmpdir, type)
        ensure_dir_exists(dir)
        self._close_current()
        filename = '%s_pred.txt' % name
        result_file = os.path.join(dir, filename)
        tar_name = os.path.join(type, filename)
        self.f = open(result_file, 'w')
        self.result_files.append((result_file, tar_name))

    def _close_current(self):
        if self.f:
            self.f.close()
            self.f = None

    def write_prediction(self, seq, interval):
        """
        Write the site information in the format expected by this test suite.

        e.g. 'M01032-16-embl, 35, 42, GGCCAAGG'
        """
        if isinstance(seq, str):
            seq = self.seq_indices[seq]
        seq_name = self.seq_names[seq]
        sequence = self.seq_dict[seq_name]
        print >> self.f, "%s, %d, %d, %s" % (
            seq_name,
            interval.start + 1,  # 1-based indexing
            interval.end,     # inclusive end-point
            sequence[interval.as_slice()].seq
        )

    def finalise(self):
        "Called when completed."
        self._close_current()
        tar = tarfile.open(self.file, "w:gz")
        for result_file, tar_name in self.result_files:
            tar.add(result_file, arcname=tar_name)
        tar.close()
        for result_file, tar_name in self.result_files:
            os.remove(result_file)
        shutil.rmtree(self.tmpdir)
        del self.tmpdir


if '__main__' == __name__:
    type = 'model_real'
    for data_set in get_data_sets_for_type(type):
        print data_set, get_fasta_file_for_data_set(type, data_set)
