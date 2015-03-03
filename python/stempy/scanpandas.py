#!/usr/bin/env python
#
# Copyright John Reid 2014
#

"""
Handle motif scans using pandas data structures.
"""

import logging
logger = logging.getLogger(__name__)

import pandas as pds
import numpy as npy


def loadoccurrences(filename='steme-pwm-scan.out'):
    """Load occurrences."""
    return pds.read_csv(
        filename,
        skiprows=1,
        header=None,
        names=(
            'motif', 'wmer', 'seq', 'pos', 'strand', 'Z', 'score', 'pvalue'))


def _renamecolumn(col):
    "Rename a column from a seq-centric scan."
    if col.startswith('# '):
        return col[2:]
    return col


def loadseqcentric(filename='scan-seq-centric.csv'):
    """Load sequence-centric statistics."""
    df = pds.read_csv(
        filename,
        header=0)
    df.columns = map(_renamecolumn, df.columns)
    return df


def loadseqinfo(filename='steme-pwm-scan.seqs'):
    """Load sequence information."""
    return pds.read_csv(
        filename,
        skiprows=1,
        header=None,
        names=('length', 'name'))


def parseseqnames(seqs):
    """Parse sequence names and add 'chrom', 'start' and 'end' columns
    to data frame."""
    colonsplit = seqs['name'].str.split(':').apply(pds.Series)
    dashsplit = colonsplit[1].str.split('-').apply(pds.Series) \
        .fillna(0).astype(npy.int)
    seqs['chrom'] = colonsplit[0].fillna("")
    seqs['start'] = dashsplit[0]
    seqs['end'] = dashsplit[1]
    return seqs


def writeseqsasbed(sequences, path_or_buf):
    """Write the sequences as a BED file."""
    sequences[sequences.chrom != ""][['chrom', 'start', 'end']].to_csv(
        path_or_buf, sep='\t', header=False, index=False)


def loadpairs(filename='spacings.out'):
    return pds.read_table(
        filename,
        sep=r'\s+',
        header=False,
        names=('primary', 'secondary', 'samestrand', 'updown', 'dist',
               'seq', 'pos', 'strand'))
