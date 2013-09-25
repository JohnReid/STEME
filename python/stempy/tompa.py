#
# Copyright John Reid 2009, 2010
#

"""
Package to handle Tompa benchmark data set.

U{http://bio.cs.washington.edu/assessment/}
"""

import os
import glob
from stempy import _read_fasta

_data_dir = '/home/john/Data/Tompa-data-set'


def get_data_set_types():
    "@return: The types of data sets."
    return ('Generic', 'Real', 'MChain')


def get_type_of_data_set(data_set):
    """
    @return: The type of the data set.
    """
    if data_set.endswith('g'):
        return 'Generic'
    if data_set.endswith('r'):
        return 'Real'
    if data_set.endswith('m'):
        return 'MChain'
    raise RuntimeError('Could not get type of data set (%s)' % data_set)


def get_directory_for_type(t):
    "@return: The directory for the data set type, t"
    return os.path.join(_data_dir, t)


def get_fasta_file_for_data_set(dataset):
    "@return: The FASTA filename for the data set."
    data_set, t = dataset
    return os.path.join(get_directory_for_type(t), '%s.fasta' % data_set)


def get_data_set_from_fasta(fasta):
    "Parse the FASTA filename and return the data set and type."
    dir, basename = os.path.split(fasta)
    data_set, ext = os.path.splitext(basename)
    assert ext == '.fasta'
    t = os.path.basename(dir)
    assert t == get_type_of_data_set(data_set)
    return data_set, t


def get_data_sets_for_type(t):
    "@return: A generator that yields the (data set, type) tuples for all data sets of the type"
    import glob
    for fasta in glob.glob(os.path.join(get_directory_for_type(t), '*.fasta')):
        yield get_data_set_from_fasta(fasta)


class ResultsFormatter(object):

    "Formats binding site predictions into a format suitable for submitting to the web service."

    def __init__(self, method_name):
        "Initialise."
        self.f = open('Tompa-%s-results.txt' % method_name, 'w')

    def start_data_set(self, dataset, fasta):
        "Called when starting to format results of particular dataset."
        # need sequences for their lengths later on...
        self.seq_dict, self.seq_indices, self.seq_names = _read_fasta(fasta)
        data_set, _t = dataset
        print >> self.f, ">dataset\n%s\n>instances" % data_set

    def write_prediction(self, seq, interval):
        "Writes prediction to results."
        seq = seq.strip()
        if isinstance(seq, str):
            seq = self.seq_indices[seq]
        sequence = self.seq_dict[self.seq_names[seq]]
        print >> self.f, "%d,%d,%s" % (
            seq, interval.start - len(sequence), sequence[interval.as_slice()].seq)

    def finalise(self):
        "Called when completed."
        self.f.close()
        del self.f


if '__main__' == __name__:
    for data_set, t in get_data_sets_for_type('Real'):
        print data_set, get_fasta_file_for_data_set(data_set)
