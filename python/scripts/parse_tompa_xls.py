#
# Copyright John Reid 2010
#

"""
Parse xls output from Tompa assessment website.
"""

import sys
import numpy as np
import pylab as pl
from elementtidy import TidyHTMLTreeBuilder
from infpy.roc import RocCalculator
from cookbook.pylab_utils import set_rcParams_for_latex


def child_texts(node):
    "@return: The text of the children of the node."
    return [entry.text for entry in node.getchildren() if entry.text]


def parse_xls_file(xls_file):
    "Parse a xls file."
    tree = TidyHTMLTreeBuilder.parse(xls_file)
    root = tree.getroot()
    table = root.find(
        '{http://www.w3.org/1999/xhtml}body/{http://www.w3.org/1999/xhtml}table')
    col_names = None
    data_sets = []
    matrix = []
    for row in table.findall('{http://www.w3.org/1999/xhtml}tr'):
        row_texts = child_texts(row)
        if None == col_names:
            col_names = row_texts[1:]
        elif len(row_texts) == len(col_names) + 1:
            data_sets.append(row_texts[0])
            matrix.append(map(float, row_texts[1:]))
    matrix = np.array(matrix)
    return data_sets, col_names, matrix


def rocs_for_row(row):
    "@return: A ROC for each of the nucleotide and site level statistics."
    nucleotide_roc = RocCalculator(row[0], row[1], row[3], row[2])
    site_roc = RocCalculator(row[4], row[5], 0, row[6])
    return nucleotide_roc, site_roc


def statistics_for_row(row):
    "@return: A list of statistics and their names."
    nroc, sroc = rocs_for_row(row)
    statistics = []
    names = []

    # Nucleotide sensitivity
    statistics.append(nroc.sensitivity())
    names.append('nSn')

    # Nucleotide positive predictive value
    statistics.append(nroc.positive_predictive_value())
    names.append('nPPV')

    # Nucleotide performance coefficient
    statistics.append(nroc.performance_coefficient())
    names.append('nPC')

    # Nucleotide correlation coefficient
    statistics.append(nroc.correlation_coefficient())
    names.append('nCC')

    # Site sensitivity
    statistics.append(sroc.sensitivity())
    names.append('sSn')

    # Site positive predictive value
    statistics.append(sroc.positive_predictive_value())
    names.append('sPPV')

    # Site average performance
    statistics.append(sroc.average_performance())
    names.append('sASP')

    return statistics, names


def make_bar_chart(meme_statistics, stem_statistics, stat_names):
    "Make a bar chart comparing MEME to STEM."
    xlocations = np.arange(len(stat_names)) + 0.5
    width = 0.25
    pl.bar(xlocations - width / 2, meme_statistics,
           width=width, label='MEME', color='green')
    pl.bar(xlocations + width / 2, stem_statistics,
           width=width, label='STEM', color='blue')
    pl.xticks(xlocations + width / 2, stat_names)
    pl.xlim(0, xlocations[-1] + width * 2)
    pl.gca().get_xaxis().tick_bottom()
    pl.gca().get_yaxis().tick_left()

#
# Get arguments
#
if len(sys.argv) != 3:
    raise RuntimeError(
        'USAGE: %s <MEME xls filename> <STEME xls filename>' % sys.argv[0])
meme_xls_file = sys.argv[1]
stem_xls_file = sys.argv[2]

#
# Parse xls files
#
data_sets, col_names, meme_matrix = parse_xls_file(meme_xls_file)
_data_sets, _col_names, stem_matrix = parse_xls_file(stem_xls_file)
assert meme_matrix.shape == stem_matrix.shape
assert data_sets == _data_sets
assert col_names == _col_names

#
# Initialise pylab
#
set_rcParams_for_latex()

#
# Make some bar charts for single organisms and for total
#
for i in xrange(-5, 0):
    data_set = data_sets[i]
    meme_statistics, stat_names = statistics_for_row(meme_matrix[i])
    stem_statistics, _stat_names = statistics_for_row(stem_matrix[i])
    assert stat_names == _stat_names
    pl.figure()
    make_bar_chart(meme_statistics, stem_statistics, stat_names)
    pl.title(data_set)
    pl.legend()
    pl.savefig('stats-%s.eps' % data_set)
    pl.clf()
