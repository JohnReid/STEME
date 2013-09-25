#
# Copyright John Reid 2012
#


"""
Code to produce minimal MEME formatted output from STEME.
"""


import numpy as N
import sys
import stempy
from itertools import imap


class MinimalMemeOutput(object):

    """
    Produces minimal MEME formatted output (http://meme.sdsc.edu/meme/doc/meme-format.html#min_format).
    """

    def __init__(self, f):
        self.f = f
        "The output file handle."

    def initialise(self, algorithm):
        """
        Write the version
        """
        # need to add revcomp to command summary as BioPython looks for it to
        # determine if using both strands
        self.f.write("""
# If you use STEME in your research, please cite:
# 
# Nucleic Acids Res. 2011 Oct;39(18):e126. Epub 2011 Jul 23.
# STEME: efficient EM to find motifs in large data sets.
# Reid JE, Wernisch L.

# Really the version of STEME
MEME version %s

ALPHABET= ACGT

strands: + -

Background letter frequencies (from dataset with add-one prior applied):
A %.3f C %.3f G %.3f T %.3f

# Command line summary:
# %s : (revcomp)
""" % (
                     stempy.__release__,
                     algorithm.freqs_with_pseudo_counts.freq(0),
                     algorithm.freqs_with_pseudo_counts.freq(1),
                     algorithm.freqs_with_pseudo_counts.freq(2),
                     algorithm.freqs_with_pseudo_counts.freq(3),
                     ' '.join(sys.argv),
                     ))

    def found_motif(self, algorithm, motif, seconds_taken):
        self.f.write("""
MOTIF STEME-%d
letter-probability matrix: alength= 4 w= %4d nsites= %3d E= %e
%s
""" % (
                     motif.idx + 1,
                     motif.model.W, motif.num_sites, N.exp(motif.log_E_value),
                     '\n'.join([
                               ' '.join([('%7.6f' % x) for x in row])
                               for row in N.exp(motif.model.bs.pssm.log_probs.values())
                               ])
                     ))

    def finalise(self, algorithm):
        pass
