#
# Copyright John Reid 2011
#


"""
Code to produce MEME-like output from STEME.
"""


import numpy as N
import os
import sys
from itertools import imap


class MemeLikeOutput(object):

    """
    Produces MEME-like output that can be parsed by biopython, bioperl or MAST parsers.
    """

    def __init__(self, f):
        self.f = f
        "The output file handle."

    def initialise(self, algorithm):
        """
        Write the version
        """
        self.f.write("""
********************************************************************************
STEME - Motif discovery tool
********************************************************************************
MEME version 4.5.0
STEME version %d.%d.%d

For further information on how to interpret these results or to get
a copy of the MEME software please access http://meme.nbcr.net.

This file may be used as input to the MAST algorithm for searching
sequence databases for matches to groups of motifs.  MAST is available
for interactive use and downloading at http://meme.nbcr.net.
********************************************************************************

""" % algorithm.version())

        #
        # Write the reference info
        #
        self.f.write("""
********************************************************************************
REFERENCE
********************************************************************************
If you use this program in your research, please cite:

Nucleic Acids Res. 2011 Oct;39(18):e126. Epub 2011 Jul 23.
STEME: efficient EM to find motifs in large data sets.
Reid JE, Wernisch L.
********************************************************************************
""")

        #
        # Write the sequence info if we don't have too many
        #
        if len(algorithm.input_sequences.seqs) < algorithm.options.max_seqs_to_write:
            self.f.write("""
********************************************************************************
TRAINING SET
********************************************************************************
DATAFILE= %s
ALPHABET= ACGT
Sequence name            Weight Length  Sequence name            Weight Length  
-------------            ------ ------  -------------            ------ ------
%s
********************************************************************************

""" % (
                algorithm.fasta,
                '\n'.join("%-24s 1.0000 %6d" % (_id, algorithm.input_sequences.seqs.string_length(i))
                          for i, _id in enumerate(algorithm.input_sequences.ids))
            ))

        #
        # Write the command line summary
        #
        # need to add revcomp to command summary as BioPython looks for it to
        # determine if using both strands
        self.f.write("""
********************************************************************************
COMMAND LINE SUMMARY
********************************************************************************
This information can also be useful in the event you wish to report a
problem with the STEME software.

command: %s : (revcomp) 

model:  mod=           anr    nmotifs= %9d    evt=           inf
object function=  E-value of product of p-values
width:  minw= %12d    maxw= %12d    minic=        0.00
nsites: minsites= %8d    maxsites= %8d    wnsites=       %.1f
theta:  prob=            1    spmap=         uni    spfuzz= %10.1f
global: substring=     yes    branching=      no    wbranch=        no
em:     prior=   dirichlet    b=            0.01    maxiter= %9d
        distance=    %.0e
data:   n= %15d    N= %15d
strands: + -
sample: seed=            0    seqfrac=         1
Letter frequencies in dataset:
A %.3f C %.3f G %.3f T %.3f 
Background letter frequencies (from dataset with add-one prior applied):
A %.3f C %.3f G %.3f T %.3f 
********************************************************************************

""" % (
                     ' '.join(sys.argv),
                     algorithm.options.num_motifs,
                     algorithm.options.min_w,
                     algorithm.options.max_w,
                     algorithm.options.min_num_sites,
                     algorithm.options.max_num_sites,
                     algorithm.options.wnsites,
                     algorithm.options.starts_seed_pseudo_counts,
                     algorithm.options.max_iters,
                     algorithm.options.convergence_distance,
                     algorithm.input_sequences.data.N,
                     len(algorithm.input_sequences.seqs),
                     algorithm.freqs.freq(0),
                     algorithm.freqs.freq(1),
                     algorithm.freqs.freq(2),
                     algorithm.freqs.freq(3),
                     algorithm.freqs_with_pseudo_counts.freq(0),
                     algorithm.freqs_with_pseudo_counts.freq(1),
                     algorithm.freqs_with_pseudo_counts.freq(2),
                     algorithm.freqs_with_pseudo_counts.freq(3),
                     ))

    def write_motif(self, algorithm, motif):
        """
        Write the motif
        """
        rounded_pspm = (
            N.exp(motif.model.bs.pssm.log_probs.values()) * 10.).round().astype(int)

        def char(x):
            if 0 == x:
                return ':'
            if 10 == x:
                return 'a'
            return str(x)

        def pssm_line(b):
            return ''.join(imap(char, rounded_pspm[:, b]))

        self.f.write("""
********************************************************************************
MOTIF %2d    width = %4d   sites = %3d   llr = %f   E-value = %e
********************************************************************************
--------------------------------------------------------------------------------
    Motif %d Description
--------------------------------------------------------------------------------
Simplified        A  %s
pos.-specific     C  %s
probability       G  %s
matrix            T  %s

--------------------------------------------------------------------------------
""" % (motif.idx + 1, motif.model.W, motif.num_sites, motif.LLR, N.exp(motif.log_E_value), motif.idx + 1,
                    pssm_line(0), pssm_line(1), pssm_line(2), pssm_line(3)))

    def write_relative_entropy(self, algorithm, motif):
        """
        Write the relative entropy - NOT IMPLEMENTED
        """
        raise NotImplementedError()
        self.f.write("""
         bits    2.1         
                 1.9 * ** * *
                 1.7 * ** * *
                 1.4 * ** * *
Relative         1.2 **** * *
Entropy          1.0 **** ***
(12.5 bits)      0.8 **** ***
                 0.6 **** ***
                 0.4 ********
                 0.2 ********
                 0.0 --------
""")

    def write_multilevel_consensus(self, algorithm, motif):
        """
        Write the multilevel consensus sequence - NOT IMPLEMENTED
        """
        raise NotImplementedError()
        self.f.write("""
Multilevel           CAGCCCTG
consensus             T  A A 
sequence                 T                             
""")

    def write_sites(self, algorithm, motif):
        """
        Write something like...
        
        FR0000854                    -    254  1.81e-05 GAGTCCACAC CAGCCCTG CCCCAGGCTC
        """
        if len(algorithm.input_sequences.seqs) < algorithm.options.max_sites_to_write:
            self.f.write("""
--------------------------------------------------------------------------------
    Motif %d sites sorted by position p-value
--------------------------------------------------------------------------------
Sequence name            Strand  Start   P-value              Site
-------------            ------  ----- ---------            --------
%s
--------------------------------------------------------------------------------
""" % (
                motif.idx + 1,
                '\n'.join([self._output_for_prediction(algorithm, prediction)
                          for prediction in motif.predictions])
            ))

    def _output_for_prediction(self, algorithm, prediction, site_flanks=10):
        """
        Sequence name            Strand  Start   P-value              Site
        -------------            ------  ----- ---------            --------
        FR0000854                    -    254  1.81e-05 GAGTCCACAC CAGCCCTG CCCCAGGCTC
        """
        from . import reverse_complement
        seq_index = algorithm.input_sequences.index_for_seq_id(prediction.seq)
        flank_start = max(prediction.interval.start - site_flanks, 0)
        flank_end = min(prediction.interval.end + site_flanks,
                        algorithm.input_sequences.data.seq_length(seq_index))
        left_flank = algorithm.input_sequences.data.subsequence(
            seq_index, flank_start, prediction.interval.start - flank_start)
        site = algorithm.input_sequences.data.subsequence(
            seq_index, prediction.interval.start, prediction.interval.end - prediction.interval.start)
        right_flank = algorithm.input_sequences.data.subsequence(
            seq_index, prediction.interval.end, flank_end - prediction.interval.end)
        if prediction.rev_comp:
            left_flank, right_flank = reverse_complement(
                right_flank), reverse_complement(left_flank)
            site = reverse_complement(site)
        if not left_flank:
            left_flank = '.'
        if not right_flank:
            right_flank = '.'
        return "%-24s %5s %6d %9.1e %*s %s %-*s" % (
            # only print first part of sequence name, otherwise parsing can get
            # confused later
            prediction.seq.split()[0],
            prediction.rev_comp and '-' or '+',
            prediction.interval.start,
            prediction.p_value,
            site_flanks, left_flank,
            site,
            site_flanks, right_flank,
        )

    def write_block_diagrams(self, algorithm, motif):
        """
        Write block diagrams - NOT IMPLEMENTED
        """
        raise NotImplementedError()
        self.f.write("""
--------------------------------------------------------------------------------
    Motif %d block diagrams
--------------------------------------------------------------------------------
SEQUENCE NAME            POSITION P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
FR0000854                         0.00012  156_[-1]_89_[-1]_18_[-1]_416_[-1]_
                                           71_[+1]_65
FR0000815                         1.8e-05  132_[-1]_4_[+1]_112_[+1]_561
FR0000878                          0.0001  364_[+1]_[-1]_23
FR0000873                         5.1e-05  305_[+1]_356
FR0000868                         0.00012  273_[+1]_20_[+1]_200_[+1]_678
FR0000822                         5.1e-05  [-1]_73
FR0000068                         5.1e-05  289_[-1]_79
FR0000841                          0.0001  88_[+1]_138_[-1]_299
FR0000804                          0.0002  150_[-1]_188_[-1]_514
FR0000534                          0.0001  305_[+1]_179
FR0000781                         0.00012  248_[+1]_19
FR0000861                          0.0002  161_[-1]_420
--------------------------------------------------------------------------------
""" % (motif.idx + 1,))

    def write_blocks(self, algorithm, motif):
        """
        Write blocks format - NOT IMPLEMENTED
        """
        raise NotImplementedError()
        self.f.write("""
--------------------------------------------------------------------------------
    Motif %d in BLOCKS format
--------------------------------------------------------------------------------
BL   MOTIF %d width=8 seqs=23
FR0000854                (  254) CAGCCCTG  1 
FR0000854                (  157) CAGCCCTG  1 
FR0000815                (  265) CAGCCCTG  1 
FR0000815                (  145) CAGCCCTG  1 
FR0000815                (  133) CAGCCCTG  1 
FR0000878                (  365) CAGCACTG  1 
FR0000873                (  306) CAGCACTG  1 
FR0000868                (  302) CAGCTCTG  1 
FR0000868                (  274) CAGCACTG  1 
FR0000822                (    1) CAGCTCTG  1 
FR0000068                (  290) CAGCTCTG  1 
FR0000854                (  783) CAGCCCAG  1 
FR0000878                (  373) CAGCACAG  1 
FR0000854                (  280) CAGCACAG  1 
FR0000841                (  235) CAGCACAG  1 
FR0000841                (   89) CAGCACAG  1 
FR0000804                (  151) CAGCTCAG  1 
FR0000534                (  306) CAGCTCAG  1 
FR0000868                (  510) CTGCCCTG  1 
FR0000854                (  704) CTGCCCTG  1 
FR0000781                (  249) CTGCCCTG  1 
FR0000861                (  162) CTGCTCAG  1 
FR0000804                (  347) CTGCTCAG  1 
//

--------------------------------------------------------------------------------
""" % (motif.idx + 1, motif.idx + 1))

    def write_pssm(self, algorithm, motif):
        """
        Write position-specific-scoring-matrix
        """
        int_log_odds = (
            motif.model.bs.pssm.log_probs.values() * 100).round().astype(int)
        self.f.write("""
--------------------------------------------------------------------------------
    Motif %d position-specific scoring matrix
--------------------------------------------------------------------------------
log-odds matrix: alength= 4 w= %2d n= %5d bayes= %8.5f E= %.1e 
%s
--------------------------------------------------------------------------------
""" % (
                motif.idx + 1, motif.model.W, motif.model.num_W_mers,
                N.log2((1 - motif.model.lambda_) / motif.model.lambda_), N.exp(
                    motif.log_E_value),
                '\n'.join(''.join('%6d' % l for l in int_log_odds_col)
                          for int_log_odds_col in int_log_odds)
            )
        )

    def write_pspm(self, algorithm, motif):
        """
        Write position-specific-probability-matrix
        """
        pspm = N.exp(motif.model.bs.pssm.log_probs.values())
        self.f.write("""
--------------------------------------------------------------------------------
    Motif %d position-specific probability matrix
--------------------------------------------------------------------------------
letter-probability matrix: alength= 4 w= %2d nsites= %3d E= %.1e 
%s
--------------------------------------------------------------------------------
""" % (
    motif.idx + 1, motif.model.W, motif.num_sites, N.exp(motif.log_E_value),
    '\n'.join(''.join('%10f' % l for l in pspm_col) for pspm_col in pspm)
))

    def write_regex(self, algorithm, motif):
        """
        Not implemented
        """
        raise NotImplementedError()
        self.f.write("""
--------------------------------------------------------------------------------
    Motif %d regular expression
--------------------------------------------------------------------------------
C[AT]GC[CAT]C[TA]G
--------------------------------------------------------------------------------
""" % (motif.idx + 1,))

    def found_motif(self, algorithm, motif, seconds_taken):
        self.write_motif(algorithm, motif)
        # self.write_relative_entropy(algorithm, motif) # not implemented yet
        # self.write_multilevel_consensus(algorithm, motif) # not implemented
        # yet
        self.write_sites(algorithm, motif)
        # self.write_block_diagrams(algorithm, motif) # not implemented yet
        # self.write_blocks(algorithm, motif) # not implemented yet
        self.write_pssm(algorithm, motif)
        self.write_pspm(algorithm, motif)
        # self.write_regex(algorithm, motif) # not implemented yet
        self.f.write("""
Time  %f secs.

********************************************************************************
""" % seconds_taken)
        # ensure that all internal buffers associated with self.f are written
        # to disk
        self.f.flush()
        os.fsync(self.f.fileno())

    def write_combined_block_diagrams(self, algorithm, motifs):
        """
        Write combined block diagrams - NOT IMPLEMENTED
        """
        raise NotImplementedError()
        self.f.write("""
--------------------------------------------------------------------------------
    Combined block diagrams: non-overlapping sites with p-value < 0.0001
--------------------------------------------------------------------------------
SEQUENCE NAME            COMBINED P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
FR0000068                        1.82e-01  376
FR0000253                        1.56e-01  170
FR0000534                        2.32e-01  492
FR0000679                        1.75e-01  360
--------------------------------------------------------------------------------
""")

    def finalise(self, algorithm):
        """
        Write an empty summary block. Need to implement this as pycogent does
        parse this block.
        """
        _sysname, nodename, _release, _version, _machine = os.uname()
        self.f.write("""

********************************************************************************
SUMMARY OF MOTIFS
********************************************************************************

--------------------------------------------------------------------------------
    Combined block diagrams: non-overlapping sites with p-value < 0.0001
--------------------------------------------------------------------------------
SEQUENCE NAME            COMBINED P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
--------------------------------------------------------------------------------

********************************************************************************


********************************************************************************
Stopped because nmotifs = %d reached.
********************************************************************************

CPU: %s

********************************************************************************
""" % (algorithm.options.num_motifs, nodename))
        self.f.close()
