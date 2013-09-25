#
# Copyright John Reid 2009, 2010, 2011, 2012, 2013
#

from __future__ import with_statement

import pkg_resources

__doc__ = pkg_resources.resource_string(__name__, "README")
__license__ = pkg_resources.resource_string(__name__, "LICENSE")
__release__, __svn_revision__ = pkg_resources.resource_string(
    __name__, "VERSION").strip().split('-')
__major_version__, __minor_version__, __release_version__ = map(
    int, __release__.split('.'))
__version__ = '%d.%d' % (__major_version__, __minor_version__)


def version_string():
    """Return the release and svn revision as a string."""
    return '%s %s' % (__release__, __svn_revision__)


from ._stempy_build import *
from ._stempy_build import _dummy_fn, _debug, _python_debug_build, _has_google_profiler
if _has_google_profiler:
    from ._stempy_build import __google_profiler_start, __google_profiler_stop
from ._meme_like_output import MemeLikeOutput
from ._minimal_meme_output import MinimalMemeOutput

from cookbook.named_tuple import namedtuple
from cookbook.timer import SimpleTimer, Timer
from cookbook.interval import Interval
from cookbook.pylab_utils import pylab_ioff, pylab_context_ioff

import glob
import os
import math
import numpy as N
import itertools
import time
import subprocess
import pylab as P
import warnings
import cPickle
from optparse import OptionParser, OptionGroup
from itertools import imap
from collections import defaultdict

logger = logging.getLogger(__name__)
try:
    logging.captureWarnings(True)
except AttributeError:
    warnings.warn(
        'Could not ask logging system to capture warnings. Old logging module?')


def add_options(option_parser):
    """
    Add options for STEME to the option parser.
    """

    #
    # Multiple motif options
    #
    multi_motif_options = OptionGroup(
        option_parser,
        "Multiple motifs",
        "Control how to find more than one motif."
    )
    multi_motif_options.add_option(
        "--num-motifs",
        dest="num_motifs",
        default=1,
        type='int',
        help="Number of motifs to look for."
    )
    multi_motif_options.add_option(
        "--prediction-Z-threshold",
        default=.3,
        type='float',
        help="The threshold on Z used to erase instances of motifs. The lower this is, the more instances will be erased."
    )
    option_parser.add_option_group(multi_motif_options)

    #
    # Output options
    #
    output_options = OptionGroup(
        option_parser,
        "Output",
        "Control the output location, format, writing logos, etc..."
    )
    output_options.add_option(
        "--output-dir",
        dest="output_dir",
        default='output',
        help="Output directory."
    )
    output_options.add_option(
        "--meme-like-output",
        default="meme.txt",
        help="Produce MEME-like output so that it can be parsed by downstream tools."
    )
    output_options.add_option(
        "--html-output",
        default="STEME.html",
        help="Produce HTML output."
    )
    output_options.add_option(
        "--print-sites",
        default=False,
        action="store_true",
        help="Write a file containing the sites that were used to make the motif."
    )
    output_options.add_option(
        "--dont-write-logos",
        default=False,
        action="store_true",
        dest="dont_write_logos",
        help="Don't write logos for motifs."
    )
    output_options.add_option(
        "--write-em-logos",
        default=False,
        action="store_true",
        dest="write_em_logos",
        help="Write logos for motifs during EM algorithm."
    )
    output_options.add_option(
        "--write-em-stats",
        default=False,
        action="store_true",
        dest="write_em_stats",
        help="Write statistics for EM algorithm."
    )
    output_options.add_option(
        "--tomtom",
        default=[],
        action="append",
        help="Run TOMTOM tool from the the MEME suite on the motifs using the specified motif databases."
    )
    output_options.add_option(
        "--max-seqs-to_write",
        default=1000,
        type="int",
        help="Maximum number of sequences to write information about in output."
    )
    output_options.add_option(
        "--max-sites-to_write",
        default=1000,
        type="int",
        help="Maximum number of sites to write information about in output."
    )
    option_parser.add_option_group(output_options)

    #
    # Background options
    #
    background_options = OptionGroup(
        option_parser,
        "Background model",
        "Control the background model."
    )
    background_options.add_option(
        "--bg-model-order",
        dest="bg_model_order",
        default=2,
        type='int',
        help="Order of the background Markov model."
    )
    background_options.add_option(
        "--bg-fasta-file",
        default=None,
        help="If specified, STEME builds its background model from the sequences in this file rather than from the input sequences."
    )
    background_options.add_option(
        "--back-dist-prior",
        dest="back_dist_prior",
        default=1.,
        type='float',
        help="Pseudo-counts for Markov background model."
    )
    option_parser.add_option_group(background_options)

    #
    # Start finding options
    #
    start_finding_options = OptionGroup(
        option_parser,
        "Start finding",
        "Control how starts are found."
    )
    start_finding_options.add_option(
        "--max-start-finding-time",
        default=0.,
        type='float',
        help="How many seconds to dedicate to finding starts for each motif. If not given, STEME will look at each possible start (can be slow)."
    )
    start_finding_options.add_option(
        "--min-sites",
        dest="min_num_sites",
        default=None,
        type='int',
        help="Minimum number of sites. Defaults to # sequences / 10."
    )
    start_finding_options.add_option(
        "--max-sites",
        dest="max_num_sites",
        default=None,
        type='int',
        help="Maximum number of sites. Defaults to 50% more than # sequences."
    )
    start_finding_options.add_option(
        "-w",
        dest="width",
        default=[],
        type='int',
        action='append',
        help="If specified, search for motifs of this width (can specify more than one)."
    )
    start_finding_options.add_option(
        "--minw",
        dest="min_w",
        default=6,
        type='int',
        help="Minimum width of motif to find."
    )
    start_finding_options.add_option(
        "--maxw",
        dest="max_w",
        default=14,
        type='int',
        help="Maximum width of motif to find."
    )
    start_finding_options.add_option(
        "--starts-per-motif",
        default=4,
        type='int',
        help="Number of starts to find per motif."
    )
    start_finding_options.add_option(
        "--use-seed",
        dest="use_seed",
        default=None,
        help="If specified, only use this seed as a start."
    )
    start_finding_options.add_option(  # same as -spfuzz MEME option
        "--starts-seed-pseudo-counts",
        dest="starts_seed_pseudo_counts",
        default=.5,
        type='float',
        help="Pseudo counts with which to smooth possible starting seeds."
    )
    start_finding_options.add_option(
        "--starts-speed-up",
        default=0,
        type='int',
        help="Speed up the start finding by ignoring so many potential starting points."
    )
    start_finding_options.add_option(
        "--candidate-starts-factor",
        default=N.sqrt(2.),
        type='float',
        help="The factor for the geometric progression that determines which numbers of sites to try when start finding."
    )
    option_parser.add_option_group(start_finding_options)

    #
    # EM options
    #
    EM_options = OptionGroup(
        option_parser,
        "EM",
        "Control the behaviour of the Expectation Maximization algorithm."
    )
    EM_options.add_option(
        "--max-iterations",
        dest="max_iters",
        default=1000,
        type='int',
        help="Maximum number of iterations for EM algorithm."
    )
    EM_options.add_option(
        "--dont-discretize",
        default=False,
        action="store_true",
        dest="dont_discretize",
        help="Don't run discretisation after EM."
    )
    EM_options.add_option(
        "--convergence_distance",
        dest="convergence_distance",
        default=1e-5,
        type='float',
        help="Threshold between successive iterations at which to stop EM."
    )
    EM_options.add_option(
        "--wnsites",
        dest="wnsites",
        default=.8,
        type='float',
        help="Weight on number of sites. Used when updating lambda in EM algorithm."
    )
    EM_options.add_option(  # same as -b MEME option
        "--em-seed-pseudo-counts",
        dest="em_seed_pseudo_counts",
        default=.01,
        type='float',
        help="Pseudo counts for motif model in EM algorithm."
    )
    EM_options.add_option(
        "--epsilon",
        dest="epsilon",
        default=0.4,
        type='float',
        help="Allowed error in motif probabilities for EM algorithm."
    )
    option_parser.add_option_group(EM_options)

    #
    # Internal options
    #
    internal_options = OptionGroup(
        option_parser,
        "Internal",
        "Not normally used."
    )
    internal_options.add_option(
        "--cache-index",
        action='store_true',
        help="Save the index to disk so it does not need to be rebuilt next time."
    )
    internal_options.add_option(
        "-Q",
        dest="pvalue_table_lattice_size",
        default=0.,
        type='int',
        help="Size of the lattice used when calculating p-values. Default is 2 * max sites."
    )
    internal_options.add_option(
        "--lambda",
        dest="lambda_",
        default=0.,
        type='float',
        help="Likelihood of a binding site in the model. Set to reasonable value by default."
    )
    internal_options.add_option(
        "--alphabet-size",
        dest="alphabet_size",
        default=4,
        type='int',
        help="Number of characters in the alphabet."
    )
    internal_options.add_option(
        "--write-IC",
        default="",
        help="Write the starts' IC values to the given file."
    )
    internal_options.add_option(
        "--google-profile",
        dest="google_profile",
        action="store_true",
        default=False,
        help="Profile with the google profiler."
    )
    internal_options.add_option(
        "--store-start-results",
        action="store_true",
        default=False,
        help="Retain the start results (used for testing)."
    )
    option_parser.add_option_group(internal_options)


def get_default_options():
    "@return: The default options."
    parser = OptionParser()
    add_options(parser)
    return parser.get_default_values()


def log_option_list(name, option_list, options):
    logger.info('%s:' % name)
    for option in option_list:
        if option.dest:
            logger.info(
                '%32s: %-32s * %s',
                option.dest, str(getattr(options, option.dest)), option.help
            )
    logger.info('')


def log_option_groups(option_parser, options):
    log_option_list("Options", option_parser.option_list, options)
    for option_group in option_parser.option_groups:
        log_option_list(option_group.title, option_group.option_list, options)


def parse_options(additional_options_adder=None, option_parser=None):
    "Parse command line options."
    if None == option_parser:
        option_parser = OptionParser()
    if additional_options_adder:
        additional_options_adder(option_parser)
    options, args = option_parser.parse_args()
    check_google_profiling(options)
    log_option_groups(option_parser, options)
    return options, args


def get_markov_model_create_fn(options):
    """Return a function that will create a Markov model of the order specified in the options.
    """
    return globals()['create_markov_model_order_from_index_%d' % options.bg_model_order]


def get_bg_model_from_markov_model_fn(options):
    """Return a function that will create a background model from a Markov model of the order specified in the options.
    """
    return globals()['create_bg_model_from_Markov_model_%d' % options.bg_model_order]


def complement(b):
    if 'A' == b:
        return 'T'
    if 'C' == b:
        return 'G'
    if 'G' == b:
        return 'C'
    if 'T' == b:
        return 'A'
    if 'N' == b:
        return 'N'
    raise RuntimeError('Unknown base: %s' % b)


def reverse_complement(s):
    return ''.join(imap(complement, s[::-1]))


def ensure_dir_exists(d):
    if not os.path.exists(d):
        os.makedirs(d)


def basename_wo_ext(filename):
    "@return: The basename of the file without any extension."
    return os.path.splitext(os.path.basename(filename))[0]


def empty_dir(d):
    "Empty a directory. Use with caution!"
    for f in glob.glob(os.path.join(d, '*')):
        os.remove(f)


def check_output(*popenargs, **kwargs):
    r"""Run command with arguments and return its output as a byte string.

    If the exit code was non-zero it raises a CalledProcessError.  The
    CalledProcessError object will have the return code in the returncode
    attribute and output in the output attribute.

    The arguments are the same as for the Popen constructor.  Example:

    >>> check_output(["ls", "-l", "/dev/null"])
    'crw-rw-rw- 1 root root 1, 3 Oct 18  2007 /dev/null\n'

    The stdout argument is not allowed as it is used internally.
    To capture standard error in the result, use stderr=STDOUT.

    >>> check_output(["/bin/sh", "-c",
    ...               "ls -l non_existent_file ; exit 0"],
    ...              stderr=STDOUT)
    'ls: non_existent_file: No such file or directory\n'
    """
    if 'stdout' in kwargs:
        raise ValueError('stdout argument not allowed, it will be overridden.')
    process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        raise subprocess.CalledProcessError(retcode, cmd, output=output)
    return output


def log_memory_usage():
    "Log the memory usage of the current process."
    logger.info(check_output(['ps', 'v', '-p', str(os.getpid())]))


def logo(dist, tag, d, make_png=False, make_eps=True, write_title=True, show_fineprint=True):
    "Generate a logo with the given tag in the given directory."
    import weblogolib as W
    import corebio.seq
    data = W.LogoData.from_counts(corebio.seq.unambiguous_dna_alphabet, dist)
    scale = 5.4 * 4
    options = W.LogoOptions(
        logo_title=write_title and tag or None,
        stack_width=scale,
        stack_aspect_ratio=5,
        color_scheme=W.colorscheme.nucleotide,
        show_xaxis=False,
        show_yaxis=False,
        show_fineprint=show_fineprint,
        fineprint='found by STEME %s' % __version__,
    )
    format_ = W.LogoFormat(data, options)
    ensure_dir_exists(d)
    filename = 'logo-%s' % tag
    if make_png:
        W.png_formatter(data,
                        format_, open(os.path.join(d, '%s.png' % filename), 'w'))
    if make_eps:
        W.eps_formatter(data,
                        format_, open(os.path.join(d, '%s.eps' % filename), 'w'))


def consensus_from_pssm(pssm, bg_freqs=None):
    "@return: The consensus sequence for the PSSM."
    if None != bg_freqs:
        pssm = pssm / bg_freqs
    indices = pssm.argmax(axis=1)
    bases = ['A', 'C', 'G', 'T']
    return ''.join(bases[i] for i in indices)


def hamming_distance(seq1, seq2):
    "@return: The Hamming distance between the 2 sequences."
    if len(seq1) != len(seq2):
        raise ValueError('Sequences do not have the same length.')
    result = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 != c2:
            result += 1
    return result


def index_for_base(b):
    if 'a' == b or 'A' == b:
        return 0
    if 'c' == b or 'C' == b:
        return 1
    if 'g' == b or 'G' == b:
        return 2
    if 't' == b or 'T' == b:
        return 3
    return 4


def occurrences_from_index(index):
    """
    Return an array holding the 0-order occurrence counts in the index.
    """
    def get_base_occurrences(base):
        i = Iterator(index)
        return i.goDownBase(base) and i.num_occurrences or 0
    freqs = N.array(map(get_base_occurrences, 'ACGTN'), dtype=float)
    unknown_warning_threshold = .9
    unknown_freq = freqs[4] / freqs.sum()
    if unknown_freq > unknown_warning_threshold:
        logger.warning(
            '%.1f%% of the bases in the sequences are unknown', 100. * unknown_freq)
    return freqs


def log_composition(occs):
    """Log the composition of the sequences based on the occurrence counts.
    """
    logger.info('%.1f%% of nucleotides unknown, rest distributed as:',
                float(occs[4] / occs.sum() * 100))
    logger.info('A: %4.1f%%', float(occs[0] / occs[:4].sum() * 100))
    logger.info('C: %4.1f%%', float(occs[1] / occs[:4].sum() * 100))
    logger.info('G: %4.1f%%', float(occs[2] / occs[:4].sum() * 100))
    logger.info('T: %4.1f%%', float(occs[3] / occs[:4].sum() * 100))
    logger.info('GC content: %.1f%%', float(
        occs[1:3].sum() / occs[:4].sum() * 100))


def print_theta(theta):
    for col in theta:
        print ' '.join('%.6f' % t for t in col)


def safe_x_log_x(x):
    "@return: x log(x) but replaces -inf with 0."
    l = N.log(x)
    result = x * l
    result[N.where(N.isneginf(l))] = 0.
    return result


def IC(freqs):
    "@return: The information content of the frequencies relative to uniform 0-order background."
    return safe_x_log_x(freqs).sum() + len(freqs) * N.log(4.)


def check_google_profiling(options):
    "Check that it makes sense to use google profiling if it is turned on."
    if _python_debug_build and options.google_profile:
        warnings.warn(
            'Cannot use google profiler with debug build. Turning profiler off.')
        options.google_profile = False


def turn_on_google_profiling_if_asked_for(options):
    "Start google profiling if asked to in the options."
    if options.google_profile:
        profile_file = os.path.join(options.output_dir, 'google.prof')
        logger.info('Starting google profiler. Writing to: %s', profile_file)
        __google_profiler_start(profile_file)


def turn_off_google_profiling_if_asked_for(options):
    "Stop google profiling if we have been asked to use it in the options."
    if options.google_profile:
        logger.info('Stopping google profiler.')
        __google_profiler_stop()


def print_stats(stats):
    # show stats
    print '# discarded         : %7d ( %5.1f%% )' % (
        stats.num_discarded, 100. *
        stats.num_discarded / stats.total_occurrences
    )
    print '# evaluated         : %7d ( %5.1f%% )' % (
        stats.num_evaluated, 100. *
        stats.num_evaluated / stats.total_occurrences
    )
    print '# occurrences       : %7d' % stats.total_occurrences
    print '# evaluations       : %7d' % stats.num_evaluations
    print 'Total weights       : %7.3f' % stats.evaluated_weight
    print 'Bound on discarded  : %7.3f' % stats.discarded_weight_bound
    print 'Error bound         : %7.1f%%' % (100. * stats.error_bound)


def check_log_probs(pssm_log_probs):
    "Checks log probs are in correct range."
    assert (N.exp(pssm_log_probs.values()).sum(axis=1) > 0.99).all()
    assert (N.exp(pssm_log_probs.values()).sum(axis=1) < 1.01).all()


def make_log_odds(pssm, background):
    "Calculate the log_2 odds of the PSSM relative to the background."
    if 2 != len(pssm.shape):
        raise ValueError('PSSMs must be 2-dimensional array.')
    if 1 != len(background.shape):
        raise ValueError('Background must be vector.')
    if pssm.shape[1] != len(background):
        raise ValueError('Expecting alphabet lengths to match.')
    log_odds = N.zeros((pssm.shape[0], pssm.shape[1] + 1))

    # for each column
    for pssm_col, lo_col in zip(pssm, log_odds):

        # for each base
        for b, (p, n) in enumerate(zip(pssm_col, background)):
            if 0 == n:
                lo_col[b] = 0.  # base with zero prob.
            else:
                lo_col[b] = N.log2(p) - N.log2(n)

            lo_col[-1] += n * lo_col[b]

    assert not N.isnan(lo_col).any()
    return log_odds


def scale_log_odds(log_odds, range_=100):
    "Scale the log odds and convert to integral range."
    offset = log_odds.min()
    lo_range = log_odds.max() - offset
    assert lo_range >= 0.
    if 0 == lo_range:
        raise ValueError('All log odds are the same.')
    scale = range_ / lo_range
    scaled_lo = (log_odds - offset) * scale
    # round to nearest int, not int below
    integral_scaled_lo = (scaled_lo + .5).astype(int)
    assert 0 == integral_scaled_lo.min()
    assert range_ == integral_scaled_lo.max()
    return integral_scaled_lo


def score_pssm_on_seq(pssm, seq):
    "Scores the sequence with the PSSM by adding together all the PSSM entries for the indexes in seq"
    return sum(col[b] for col, b in zip(pssm, seq))


def calc_pssm_pdf(scaled_lo, range_, background):
    "Calculate the distribution of different scores from the scaled PSSM log odds matrix."

    # set up size of arrays
    W = len(scaled_lo)
    size = W * range_ + 1

    # use 2 arrays to step through columns
    pdf_old = N.empty(size)
    pdf_new = N.zeros(size)
    pdf_new[0] = 1.

    # for each column add to distribution
    for w, scaled_lo_col in enumerate(scaled_lo):
        old_max_score = w * range_
        new_max_score = old_max_score + range_
        pdf_old, pdf_new = pdf_new, pdf_old
        pdf_new[:new_max_score + 1] = 0.
        for b, n in enumerate(background):
            score = int(scaled_lo_col[b])  # should probably be an int already
            assert 0 <= score
            assert score <= range_
            pdf_new[score:score + old_max_score + 1] += n * \
                pdf_old[:old_max_score + 1]

    # check is a proper pdf
    assert 1. - 1.e8 < pdf_new.sum()
    assert pdf_new.sum() < 1. + 1.e8

    return pdf_new


def calc_pssm_cdf(scaled_lo, range_, background):
    "Calculate the cumulative distribution of different scores from the scaled PSSM log odds matrix."
    pdf = calc_pssm_pdf(scaled_lo, range_, background)

    # compute 1-cdf from the pdf from the right to preserve right accuracy
    for i in xrange(len(pdf) - 2, -1, -1):
        pdf[i] += pdf[i + 1]

    # return the cdf
    return pdf.clip(0, 1.)


def create_log_odds(background, pssm, roundto=0, log_odds_range=100):
    """Create log odds matrix for p-values
    """
    log_odds = make_log_odds(pssm=pssm, background=background)
    scaled_log_odds = scale_log_odds(log_odds, range_=log_odds_range)
    if roundto:
        background = background.round(
            roundto)  # round background so results same as MEME (and MAST)
    pssm_cdf = calc_pssm_cdf(
        scaled_log_odds, range_=log_odds_range, background=background)
    return scaled_log_odds, pssm_cdf


InstanceInfo = namedtuple(
    'InstanceInfo', 'instance seqname seq pos wmer score pvalue')
"An instance with some extra info about the W-mer and its position."


def instance_info(instance, data, ids, W, scaled_log_odds=None, pssm_cdf=None):
    """Returns a tuple of information about the instance.
    """
    logger.debug('Getting local position: %d', instance.global_pos)
    seq, pos = data.pos_localise(instance.global_pos)
    logger.debug('Getting W-mer: %d %d', W, instance.global_pos)
    wmer = data.get_W_mer(W, instance.global_pos)
    if instance.rev_comp:
        wmer = reverse_complement(wmer)
    if None != scaled_log_odds and None != pssm_cdf:
        logger.debug('Scoring W-mer: %s', wmer)
        score = score_pssm_on_seq(scaled_log_odds, imap(index_for_base, wmer))
        pvalue = pssm_cdf[score]
    else:
        score = pvalue = None
    return InstanceInfo(
        instance=instance,
        seqname=ids[seq],
        seq=seq,
        pos=pos,
        wmer=wmer,
        score=score,
        pvalue=pvalue,
    )


def cmp_predictions(pred1, pred2):
    "Compare 2 predictions."
    seq_idx1, _seq_id1, _seq_len1, start1, site1 = pred1
    seq_idx2, _seq_id2, _seq_len2, start2, site2 = pred2
    if seq_idx1 < seq_idx2:
        return -1
    if seq_idx1 > seq_idx2:
        return 1
    if start1 < start2:
        return -1
    if start1 > start2:
        return 1
    if len(site1) < len(site2):
        return -1
    if len(site1) > len(site2):
        return 1
    return 0


def _read_fasta(fasta):
    "Read a fasta file into a dict and remember sequence indices."
    from Bio import SeqIO
    seq_dict = dict()
    seq_indices = dict()
    seq_names = dict()
    num_bases = 0
    for idx, record in zip(itertools.count(), SeqIO.parse(open(fasta, "r"), "fasta")):
        seq = record.id.strip()
        seq_dict[seq] = record
        seq_indices[seq] = idx
        seq_names[idx] = seq
        seq_length = len(record)
        num_bases += seq_length
        if seq_length > index_max_seq_length:
            warnings.warn('Sequence too long for index type: max=%d' %
                          index_max_seq_length)
    logger.info('Read %d bases from %d sequences in %s',
                num_bases, len(seq_dict), fasta)
    if seq_length > index_max_seqs:
        warnings.warn(
            'Too many sequences for index type: max=%d' % index_max_seqs)
    if num_bases > index_max_text_length:
        warnings.warn('Text too long for index type: max=%d' %
                      index_max_text_length)
    return seq_dict, seq_indices, seq_names


def read_sequences(fasta_file, options):
    """Read sequences and builds index. Will save index if options.save_index
    is set. If saved index is present will load the index instead of re-reading
    sequences and building the index.
    """
    index_file = fasta_file + '.esa'
    info_file = fasta_file + '.info'
    successful_load = False  # Did we load the index successfully?
    #if True:
    try:
        logger.info('Looking for cached index information: %s', info_file)
        num_bases, ids = cPickle.load(open(info_file))
        logger.info('Loading index of %d bases from %s.*',
                    num_bases, index_file)
        index = load_index(index_file)
        seqs = index.text()
        successful_load = True  # Did we load the index successfully?
    except:
        logger.info('Could not load cached index, reading FASTA directly')
        num_bases, seqs, ids = read_fasta(fasta_file)
        logger.info('Read %d bases from %d sequences in %s',
                    num_bases, len(seqs), fasta_file)
        logger.info('Initialising index')
        index = build_index(seqs)
        try:  # put in try-clause as not all options may have it
            if options.cache_index:
                logger.info('Creating index')
                index.visit(max_depth=options.max_w)
                logger.info('Saving index info to %s', info_file)
                cPickle.dump((num_bases, ids), open(info_file, 'w'))
                logger.info('Saving index to %s', index_file)
                save_index(index, index_file)
        except AttributeError:
            warnings.warn('Options may not have save_index attribute')
    return num_bases, seqs, ids, index


def euclidean_distance(theta1, theta2):
    "The Euclidean distance between two thetas."
    return math.sqrt(((theta1 - theta2) ** 2).sum())


def converged(last_theta, new_theta, options):
    "Is the euclidean distance between successive thetas smaller than options.convergence_distance?"
    return euclidean_distance(last_theta, new_theta) < options.convergence_distance


def initialise_random_pssm(W, alphabet_size, alpha):
    "Set up a PSSM"
    pssm_probs = N.random.dirichlet(N.ones(alphabet_size) * alpha, W)
    pssm_probs.T[:] /= pssm_probs.sum(axis=1).T
    pssm_log_probs = PssmStorage(W=W, alphabet_size=alphabet_size)
    pssm_log_probs.values()[:] = N.log(pssm_probs)
    return pssm_log_probs


def initialise_uniform_pssm(W, alphabet_size):
    "Set up a PSSM"
    pssm_probs = N.ones((W, alphabet_size)) / alphabet_size
    pssm_log_probs = PssmStorage(W=W, alphabet_size=alphabet_size)
    pssm_log_probs.values()[:] = N.log(pssm_probs)
    return pssm_log_probs


def geometric_progression(min_val, max_val, factor=2):
    "A generator for a geometric progression."
    if factor <= 1.:
        raise ValueError('factor must be > 1: %f' % factor)
    if max_val < min_val:
        raise ValueError(
            'Maximum value < minimum value : %d < %d' % (max_val, min_val))
    if min_val == int(min_val * factor):
        raise ValueError(
            'factor is too small to make a progression: %f' % factor)
    while min_val < max_val:
        yield min_val
        min_val = int(min_val * factor)
    yield max_val


def rationalise_predictions(predictions):
    "Take a sorted list of predictions and merge those that overlap."
    last_seq = last_interval = None
    for seq, interval, _rev_comp in predictions:
        if last_seq == seq and last_interval.overlap(interval):
            last_interval = last_interval.hull(interval)
        else:
            if None != last_seq:
                yield last_seq, last_interval
            last_seq, last_interval = seq, interval
    if last_seq:
        yield last_seq, last_interval


def get_Wmer(data, W, wmer):
    "@return: The W-mer, reverse complemented if necessary."
    s = data.get_W_mer(W, wmer.global_pos)
    return wmer.rev_comp and reverse_complement(s) or s


def Wmer_getter(data, W):
    return lambda wmer: get_Wmer(data, W, wmer)


def get_Wmer_info(data, W, wmer):
    "@return: info about the W-mer as a string."
    return 'Z=%.2e; pos=%9d %s; %s; %s' % (
        wmer.Z, wmer.global_pos, wmer.rev_comp and '-' or '+', data.get_W_mer(
            W, wmer.global_pos), get_Wmer(data, W, wmer)
    )


def Wmer_infoizer(data, W):
    return lambda wmer: get_Wmer_info(data, W, wmer)


def log_stats(stats):
    "Log some suffix tree statistics."
    logger.info('nodes discarded: %d', stats.nodes_discarded)
    logger.info('nodes evaluated: %d', stats.nodes_evaluated)
    logger.info('W-mers discarded: %d', stats.w_mers_discarded)
    logger.info('W-mers evaluated: %d', stats.w_mers_evaluated)


@pylab_ioff
def test_matplotlib():
    "Test matplotlib save as EPS problem."
    import pylab as P
    P.figure()
    width = .4
    x1 = [-0.4,  0.6,  1.6,  2.6,  3.6,  4.6,  5.6]
    y1 = [1.0, 0.96875, 0.77584381616639686, 0.54678532728659146,
          0.4043846450263317, 0.28361561411668612, 1.0]
    x2 = [0, 1, 2, 3, 4, 5, 6]
    y2 = [
        1.0, 0.97032301818547173, 0.77110383361936519, 0.54221305796207875,
        0.40019201884735922, 0.28326596333427007, 1.0
    ]
    P.bar(x1, y1, color='blue', width=width, label='nodes')
    P.bar(x2, y2, color='green', width=width, label='occurrences')
    P.savefig('fraction-evaluated.eps')
    P.close()


@pylab_ioff
def graph_stats(stats, tag):
    "Graph the statistics."
    import pylab as P
    # fraction evaluated plot
    P.figure()
    width = .4
    x_range = N.arange(0, len(stats.counts))
    P.bar(
        x_range - width,
        [c.node.fraction_evaluated for c in stats.counts],
        color='blue',
        width=width,
        label='nodes'
    )
    P.bar(
        x_range,
        [c.occurrence.fraction_evaluated for c in stats.counts],
        color='green',
        width=width,
        label='occurrences'
    )
    P.ylim(0, 1)
    P.xlabel('W')
    P.ylabel('Fraction evaluated')
    P.legend(loc='lower left')
    P.axes().set_xticks(range(len(stats.counts)))
    P.savefig('%s-fractions.eps' % tag)
    P.savefig('%s-fractions.png' % tag)
    P.close()

    # nodes evaluated plot
    P.figure()
    x_range = N.arange(0, len(stats.counts))
    evaluated = [c.node.evaluated for c in stats.counts]
    discarded = [c.node.discarded for c in stats.counts]
    P.bar(x_range, evaluated,       color='blue',
          label='evaluated', align='center')
    P.bar(x_range, discarded,       color='red',
          label='discarded', align='center', bottom=evaluated)
    P.xlabel('W')
    P.ylabel('Nodes')
    P.legend(loc='upper left')
    P.axes().set_xticks(range(len(stats.counts)))
    P.savefig('%s-nodes.eps' % tag)
    P.savefig('%s-nodes.png' % tag)
    P.close()

    # occurrences evaluated plot
    P.figure()
    x_range = N.arange(0, len(stats.counts))
    evaluated = [c.occurrence.evaluated for c in stats.counts]
    discarded = [c.occurrence.discarded for c in stats.counts]
    P.bar(x_range, evaluated,       color='green',
          label='evaluated', align='center')
    P.bar(x_range, discarded,       color='purple',
          label='discarded', align='center', bottom=evaluated)
    P.xlabel('W')
    P.ylabel('Occurrences')
    P.legend(loc='upper right')
    P.axes().set_xticks(range(len(stats.counts)))
    P.savefig('%s-occurrences.eps' % tag)
    P.savefig('%s-occurrences.png' % tag)
    P.close()


def log_best_w_mers(data, W, best_w_mer_finder, num_w_mers, level=logging.DEBUG, consensus=None):
    best_w_mers = list(best_w_mer_finder.best_w_mers)
    distance = '-'
    for w_mer in best_w_mers[:num_w_mers]:
        s = get_Wmer(data, W, w_mer)
        if None != consensus:
            distance = str(hamming_distance(s, consensus))
        logger.log(
            level,
            'W-mer: %s with score=%.3e on \'%s\' strand, distance=%s',
            s, w_mer.Z, w_mer.rev_comp and '-' or '+', distance
        )


EMResult = namedtuple(
    'EMResult', 'EM EM_num_sites LLs model cons_after_em cons em_duration discretise_instances coverage')
StartResult = namedtuple('StartResult', 'log_E_value start em_result')
Motif = namedtuple(
    'Motif', 'idx model num_sites LLR log_E_value start_results predictions input_stats bg_stats')
Prediction = namedtuple('Prediction', 'seq interval rev_comp p_value')
NumInstances = namedtuple('NumInstances', 'count num_bases num_seqs per_base')


def analyse_num_instances(instances, seqs):
    """
    Analyse the number of instances found in the sequences.
    """
    return NumInstances(len(instances), seqs.num_bases, len(seqs.seqs), len(instances) / float(seqs.num_bases))


def compare_instances(instances1, W1, instances2, W2):
    """
    Compare sets of instances to measure overlaps
    """
    overlap = calculate_overlap(instances1, W1, instances2, W2)
    size1 = len(instances1) * W1
    size2 = len(instances2) * W2
    return overlap, size1, size2


class StartingPointFinder(object):

    "Finds best starting points for EM."

    def __init__(self, input_sequences, bg_manager, significance, options, output_dir=None):
        "Construct."

        self.options = options
        "The options."

        self.output_dir = output_dir
        "Output directory for the start."
        if not self.output_dir:
            self.output_dir = os.path.join(self.options.output_dir, 'starts')
            ensure_dir_exists(self.output_dir)

        self.input_sequences = input_sequences
        "The input sequences."

        self.bg_manager = bg_manager

        self.significance = significance
        "The significance calculator."

        self.central_widths = self.options.width and self.options.width or list(
            geometric_progression(self.options.min_w, self.options.max_w, N.sqrt(2.)))
        "The widths for which we will find starts."
        logger.info('Looking for motifs of widths: %s',
                    ', '.join([('%d' % w) for w in self.central_widths]))

        self.start_counts = dict()
        "Counts the number of starts for each width."

        self.em_results = defaultdict(list)
        "Maps starts to EM results."

    def remove_overlapping(self, instances, W, threshold=.5):
        """
        Remove any starts that overlap the given start.
        """
        def overlaps(other):
            overlap = calculate_overlap(
                instances, W, other.best_w_mers, other.model.W)
            size = len(other.best_w_mers) * W
            if overlap > size * threshold:
                return True
            return False

        num_left = 0
        num_originally = 0
        for starts_for_width in self.starting_points:
            for starts in starts_for_width:
                num_originally += len(starts.value)
                not_overlapping = [
                    other for other in starts.value if not overlaps(other)]
                starts_for_width[starts.key][:] = not_overlapping
                num_left += len(not_overlapping)
        logger.info(
            'Left with %d/%d starts after removing overlapping.', num_left, num_originally)

    def _num_starts_for_width(self, W):
        "@return: The number of starts for the given width."
        if W not in self.start_counts:
            self.start_counts[W] = create_start_counter(
                self.input_sequences.data, W)()
        return self.start_counts[W]

    def _create_start_finder(self, W):
        "@return: start_finder"
        logger.debug('Creating start finder')
        model = self.input_sequences.create_model(
            self.bg_manager.get_bg_model(W), W)
        model.bs.seed_pseudo_counts = self.options.starts_seed_pseudo_counts
        start_finder = FindStarts(
            model.data,
            model,
            self.significance,
            self.options.min_num_sites,
            self.options.max_num_sites,
            self.options.candidate_starts_factor,
            num_to_find=self.options.starts_per_motif *
            self.options.num_motifs,
            speed_up=self.options.starts_speed_up,
        )
        start_finder.num_starts = self._num_starts_for_width(W)

        # If asked to write the information content of each start to a file,
        # then register callback to do so
        if self.options.write_IC:
            start_finder.ic_file = open(
                os.path.join(self.output_dir, self.options.write_IC.encode()), 'a')

            def write_IC_callback(instances, start, num_sites, score):
                model = start.model
                start_finder.ic_file.write(
                    'N=%d; IC=%s\n' % (
                        model.bs.pssm.num_samples,
                        ','.join(str(model.bs.pssm.columnIC(model.bg.freqs, i))
                                 for i in xrange(model.W))
                    )
                )
            start_finder.register_callback(
                FindStarts.start_examined_callback.from_callable(write_IC_callback))

# A callback that logs the W-mers
#        def Wmer_logging_callback(instances, start, num_sites, score):
# logger.info('Start finder callback: start=%s; # sites=%6d; score=%6.2f', start, num_sites, score)
#            for wmer in instances:
#                logger.info(get_Wmer_info(start.model.data, W, wmer))
#            1/0
# start_finder.register_callback(FindStarts.start_examined_callback.from_callable(Wmer_logging_callback))

        return start_finder

    def find_starts(self):
        "Find the starting points."
        if self.options.use_seed:
            logger.info('Finding starts for seed: %s', self.options.use_seed)
            self.starting_points = [
                self.starts_for_width(len(self.options.use_seed))
            ]
        else:
            logger.info(
                'Finding starting points for %d different widths: %s',
                len(self.central_widths), self.central_widths
            )
            self.starting_points = [
                self.starts_for_width(w) for w in self.central_widths
            ]

    def starts_for_width(self, W):
        "@return: The starts for the given width."
        logger.info('Examining starts of width %d', W)
        start_finder = self._create_start_finder(W)
        start_time = time.clock()
        if self.options.use_seed:
            start_finder.find_starts_for_seed(self.options.use_seed.encode())
        else:
            #
            # If we've been given an option to limit the time spent on starts then use it
            #
            if self.options.max_start_finding_time:
                self.add_time_control_callback(
                    start_finder, self.options.max_start_finding_time / len(self.central_widths))

            #
            # Find the starts
            #
            start_finder.find_starts()
            assert start_finder.start_counter == start_finder.num_starts

        logger.info(
            'Examined %d of %d starts for width %d in %f seconds',
            start_finder.starts_examined, start_finder.start_counter, W, time.clock(
            ) - start_time
        )
        for starts in start_finder.best_starts:
            for start in starts.value:
                self.log_start_info(start_finder.data, start)
        tag = os.path.join(self.output_dir, 'starts-w=%d' % W)
        graph_stats(start_finder.efficiency_statistics, tag)
        return start_finder.best_starts

    def add_time_control_callback(self, start_finder, allowed_time):
        """
        Add a callback to the start finder that controls how long it spends executing.
        This callback times how long each start takes to evaluate and changes the number that
        are ignored in order to finish on time.
        """
        total_to_examine = start_finder.num_starts
        logger.info(
            'Have been allocated %d seconds to find starts at this width, %d. Might skip some of the %d starts to achieve this.',
            allowed_time, start_finder.W, total_to_examine
        )
        start_time = time.clock()

        #
        # define the callback
        #
        def control_time_callback(instances, start, num_sites, score):
            "Our callback to control the time."
            time_taken_so_far = time.clock() - start_time
            #logger.log(log_level, 'Time taken so far: %f', time_taken_so_far)
            # if we haven't taken any time so far, then ignore...
            if time_taken_so_far > 0.:
                time_left = allowed_time - time_taken_so_far
                #logger.log(log_level, 'Time left: %f', time_left)
                if time_left <= 0:
                    # we've already used up all our allocated time
                    start_finder.speed_up = total_to_examine
                else:
                    # always under-estimate the time per start a little bit to make
                    # sure we don't skip too many
                    time_per_start = time_taken_so_far / \
                        start_finder.starts_examined
                    #logger.log(log_level, 'Time/start: %f', time_per_start)
                    how_many_to_examine = time_left / time_per_start
                    starts_left = total_to_examine - start_finder.start_counter
                    assert starts_left >= 0
                    start_finder.speed_up = int(
                        starts_left / how_many_to_examine)

        #
        # register the callback
        #
        start_finder.register_callback(
            start_finder.start_examined_callback.from_callable(control_time_callback))

    def log_start_info(self, data, start, tag=''):
        "Log some info about the start."
        self.log_start(start, self.options.max_w, tag)
        self.log_best_wmers(data, start)

    @staticmethod
    def log_start(start, max_W, tag='', level=logging.INFO):
        "Log the start."
        logger.log(
            level,
            'seed=%*s, score=%8.2f, %6d sites; %s',
            max_W, start.seed,
            start.score,
            start.num_sites,
            tag,
        )

    def log_best_wmers(self, data, start):
        w_mer_evals = list(start.best_w_mers)
        w_mer_evals.sort(reverse=True)
        w_mer_evals = w_mer_evals[:start.num_sites]
        for _eval in w_mer_evals:
            logger.debug(
                'Seed: %s; Site: %s; p(binding): %.2e; p(not binding): %.2e',
                start.seed, data.get_W_mer(
                    len(start.seed), _eval.global_pos), _eval.Z, 1. - _eval.Z
            )


class MotifFinder(object):

    """
    Finds one motif using the STEME algorithm.
    """

    def __init__(
        self,
        input_sequences,
        bg_manager,
        significance,
        ids,
        motif_idx,
        bg_freqs,
        options,
        bg_sequences=None,
        bg_bg_manager=None
    ):
        self.input_sequences = input_sequences
        "The input_sequences we are finding the motif in."

        self.bg_manager = bg_manager
        self.bg_bg_manager = bg_bg_manager

        self.significance = significance
        "Calculates the significance of motifs."

        self.ids = ids
        "The ids of the sequences we are finding the motif in."

        self.motif_idx = motif_idx
        "The index of this motif."

        self.bg_freqs = bg_freqs
        "Background frequencies (with pseudo-counts)."

        self.options = options
        "The options."

        self.bg_sequences = bg_sequences
        "The background sequences."

        # make sure output directory exists
        ensure_dir_exists(self._output_dir_for_motif())

    def _output_dir_for_motif(self):
        "@return: The output directory for this motif."
        return os.path.join(self.options.output_dir, 'motif-%02d' % self.motif_idx)

    def _output_dir_for_start(self, start):
        "@return: The directory that the output for the start is put in."
        return os.path.join(self._output_dir_for_motif(), 'start-W=%02d-%05d-%s' % (len(start.seed), start.num_sites, start.seed))

    def _run_em_from_start(self, start, callback=None, return_EM=False):
        "Run EM from a given seed (start)."
        if None == start.model:
            model = self._create_model(len(start.seed))
        else:
            model = start.model.copy()

        overall_timer = SimpleTimer()

        logger.debug(
            'Running EM from seed %s with %d sites', start.seed, start.num_sites)
        logo_tag = 'w=%02d-%s-%03d' % (
            len(start.seed), start.seed, start.num_sites)
        start_dir = self._output_dir_for_start(start)

        # seed model and set up parameters/priors
        model.bs.seed_pseudo_counts = self.options.starts_seed_pseudo_counts
        model.bs.seed(start.seed, True)
        model.prior_num_sites = start.num_sites
        model.set_lambda_for_sites(start.num_sites)
        model.bs.seed_pseudo_counts = self.options.em_seed_pseudo_counts / \
            4.  # divide by alphabet length just as MEME does.
        logger.debug(
            'Starting point for EM : IC=%3.2f, lambda=%f, pssm log sum=%f',
            IC(N.exp(model.bs.pssm.log_probs.values())),
            model.lambda_,
            model.bs.pssm.log_probs.values().sum()
        )

        # initialise start directory
        empty_dir(start_dir)
        if not self.options.dont_write_logos:
            logo(N.exp(model.bs.pssm.log_probs.values()), '%s-1-seed' %
                 logo_tag, start_dir)

        # create EM algorithm
        EM = create_EM_descender(
            self.input_sequences.data, model, self.options.epsilon, self.options.wnsites)
        EM.using_sparse_Z = False
        logger.debug(
            'EM parameters: epsilon=%f; epsilon=%f; expected_sites=%f; wnsites=%f',
            EM.epsilon, EM.epsilon, EM.expected_sites, EM.wnsites
        )
        logger.debug('EM set-up took %.3f secs', overall_timer.duration())
        overall_timer.restart()

        # run EM
        LLs = []
        timings = []
        timer = SimpleTimer()
        cumulative_duration = 0.
        theta = N.exp(model.bs.pssm.log_probs.values())
        for num_iter in xrange(self.options.max_iters):
            last_theta = theta
            timer.restart()
            expected_sites = EM.do_iteration()
            theta = N.exp(model.bs.pssm.log_probs.values())
            # print_theta(theta)
            theta_distance = euclidean_distance(last_theta, theta)
            duration = timer.duration()
            cumulative_duration += duration
            timings.append(duration)
            LLs.append(EM.LL)
            logger.debug(
                'Iteration %4d: %8.2f sites; %5.2f IC; %7d non-zero; %7d large;' +
                ' %7d normalised; sparse=%d; %.1e secs (%.1e cumulative); %.1e distance',
                num_iter, expected_sites, IC(
                    theta), EM.num_Z_non_zero, EM.num_Z_large,
                EM.num_Z_normalised, EM.using_sparse_Z, duration, cumulative_duration, theta_distance
            )
            if None != callback:
                callback(steme=self, iter=num_iter, em=EM, theta=theta,
                         theta_distance=theta_distance, duration=duration)
            if self.options.write_em_logos:
                logo(
                    N.exp(model.bs.pssm.log_probs.values()), '%s-2-iter=%03d' %
                    (logo_tag, num_iter), start_dir)

            # check convergence
            if theta_distance < self.options.convergence_distance:
                break

            # decide whether to use sparse Z for next iteration
            EM.using_sparse_Z = EM.num_Z_non_zero / \
                float(
                    self.input_sequences.data.num_occurrences(EM.W)) < 7000 / 500000.
            # EM.using_sparse_Z = num_iter % 2 # alternate using sparse Z to
            # see which is quickest
        else:
            logger.warning(
                'EM did not converge in %d iterations', num_iter + 1)

        em_duration = overall_timer.duration()
        logger.debug('EM iterations took %.3f secs', em_duration)
        overall_timer.restart()

        #logger.info("Timings: %s", ", ".join("%.1e" % t for t in timings))
        #timings_sum = sum(timings)
        #logger.info("Timings fractions: %s", ", ".join("%.1e" % (t/timings_sum) for t in timings))

        logger.debug('Lambda ratios: %s', ', '.join(('%.5f' % lr)
                     for lr in EM.lambda_ratios))

        # log result of EM
        cons_after_em = consensus_from_pssm(
            model.bs.pssm.log_probs.values(), self.bg_freqs)
        StartingPointFinder.log_start(
            start,
            self.options.max_w,
            tag='%24s -> consensus=%*s, IC=%5.1f, log10 E-value=%9.2f, %7.1f sites, LL=%.3e, lambda=%.1e' % (
                '%d EM iterations' % (num_iter + 1),
                self.options.max_w, cons_after_em,
                IC(N.exp(model.bs.pssm.log_probs.values())),
                self.significance.log_E_value(model),
                expected_sites,
                EM.LL,
                model.lambda_,
            )
        )

        if not self.options.dont_write_logos:
            logo(N.exp(model.bs.pssm.log_probs.values()), '%s-3-after-EM' %
                 logo_tag, start_dir)
        #logger.debug("Motif:\n%s" % N.exp(model.bs.pssm.log_probs.values()))
        #logger.debug("%.2f sites in model which achieves log E-value of %e", model.bs.num_samples, self.significance(model))

        # log best W-mers
    #    best_w_mer_finder = create_best_w_mer_finder(self.input_sequences.data, model, self.options.max_num_sites)
    #    best_w_mer_finder()
    #    log_best_w_mers(
    #        self.input_sequences.data,
    #        len(start.seed),
    #        best_w_mer_finder,
    #        self.options.max_num_sites,
    #        level=logging.DEBUG,
    #        consensus=consensus
    #    )

        logger.debug(
            'EM logging result took %.3f secs', overall_timer.duration())
        overall_timer.restart()

        # discretize result
        if not self.options.dont_discretize:
            _score, discretise_instances, coverage = self._discretize(
                model, start)
            logger.debug(
                'EM discretization took %.3f secs', overall_timer.duration())
            cons = consensus_from_pssm(
                model.bs.pssm.log_probs.values(), self.bg_freqs)
            overall_timer.restart()
        else:
            _score = 0.
            cons = cons_after_em
            discretise_instances = InstanceVec()

        if self.options.write_em_stats:
            graph_stats(
                EM.efficiency_statistics, os.path.join(start_dir, 'em'))

        logger.debug('EM wrap-up took %.3f secs', overall_timer.duration())
        overall_timer.restart()

        return EMResult(
            EM=return_EM and EM or None,
            EM_num_sites=expected_sites,
            LLs=LLs,
            model=model,
            cons_after_em=cons_after_em,
            cons=cons,
            em_duration=em_duration,
            discretise_instances=discretise_instances,
            coverage=coverage
        )

    def _discretize(self, model, start):
        "Look for best W-mers and build model from them. @return: The score of the model."
        # logger.debug('Discretizing model.')
        logo_tag = 'w=%02d-%s-%03d' % (
            len(start.seed), start.seed, start.num_sites)
        start_dir = self._output_dir_for_start(start)

        # functor that scores this number of sites
        def get_score_for_num_sites(num_sites_to_use):
            best_w_mer_finder.update_model(
                num_sites_to_use, use_pseudo_counts=True)
            return self.significance.log_E_value(model), num_sites_to_use

        # find the best W-mers under the model
        best_w_mer_finder = create_best_w_mer_finder(
            self.input_sequences.data, model, self.options.max_num_sites)
        best_w_mer_finder()

        # score each possible num_sites
        scores_by_num_sites = map(get_score_for_num_sites, xrange(
            self.options.min_num_sites, self.options.max_num_sites + 1))

        # sort the scores and select best
        scores_by_num_sites.sort()
        log_E_value, num_sites_to_use = scores_by_num_sites[0]

        # update model to correct motif for this num_sites - this time use
        # pseudo-counts
        best_w_mer_finder.update_model(
            num_sites_to_use, use_pseudo_counts=True)
        # as we have pseudo-counts
        assert (
            0 < N.exp(best_w_mer_finder.model.bs.pssm.log_probs.values())).all()
        # as we have pseudo-counts
        assert (
            1 > N.exp(best_w_mer_finder.model.bs.pssm.log_probs.values())).all()

        # log info
        log_10_E_value = log_E_value / N.log(10.)
        consensus = consensus_from_pssm(
            model.bs.pssm.log_probs.values(), self.bg_freqs)
        instances = InstanceVec()
        instances.extend(best_w_mer_finder.best_w_mers[:num_sites_to_use])
        sort_instances_by_position(instances)
#        logging.info('Best W-mers from start:')
#        for instance in start.best_w_mers:
#            logging.info('%6d %f %s', instance.global_pos, instance.Z, instance.rev_comp and '-' or '+')
#        logging.info('Instances:')
#        for instance in instances:
#            logging.info('%6d %f %s', instance.global_pos, instance.Z, instance.rev_comp and '-' or '+')
        overlap, start_size, _discretized_size = compare_instances(
            start.best_w_mers, start.model.W, instances, start.model.W)
        coverage = start_size and overlap / float(start_size) or 1.
        StartingPointFinder.log_start(
            start,
            self.options.max_w,
            tag='%24s -> consensus=%*s, IC=%5.1f, log10 E-value=%9.2f, coverage=%.2f' % (
                'discretised %d sites' % num_sites_to_use,
                self.options.max_w, consensus,
                IC(N.exp(best_w_mer_finder.model.bs.pssm.log_probs.values())),
                log_10_E_value,
                coverage
            )
        )

        # log best W-mers
        log_best_w_mers(
            self.input_sequences.data,
            len(start.seed),
            best_w_mer_finder,
            num_sites_to_use,
            level=logging.DEBUG,
            consensus=consensus
        )

        # write logo
        logo(N.exp(model.bs.pssm.log_probs.values()),
             '%s-4-discretised' % logo_tag, start_dir)

        return log_E_value, instances, coverage

    def _make_predictions_using_best_W_mers(self, best_model, num_sites):
        """
        Make predictions using the best so many sites.
        """
        best_w_mer_finder = create_best_w_mer_finder(
            self.input_sequences.data, best_model, num_sites)
        best_w_mer_finder()
        assert num_sites == len(best_w_mer_finder.best_w_mers)
        return best_w_mer_finder.best_w_mers

    def _make_predictions_using_instance_finder(self, best_model):
        """
        Make predictions using the best so many sites.
        """
        instance_finder = FindInstances(
            self.input_sequences.data, best_model, self.options.prediction_Z_threshold)
        instance_finder()
        return instance_finder.instances

    @pylab_ioff
    def _plot_num_sites_graph(self, start_results):
        """
        Plot information about how number of sites has changed for each start
        """
        P.figure()
        for start_result in start_results:
            P.plot([
                start_result.start.num_sites,
                start_result.em_result.EM_num_sites,
                len(start_result.em_result.discretise_instances),
            ])
        P.gca().get_xaxis().set_ticklabels(('start', 'EM', 'discretise'))
        P.gca().get_xaxis().set_ticks((0, 1, 2))
        P.ylabel('\# sites')
        num_sites_graph_filename = os.path.join(
            self._output_dir_for_motif(), 'num-sites')
        logger.info(
            'Saving number of sites graph to %s', num_sites_graph_filename)
        P.savefig('%s.eps' % num_sites_graph_filename)
        P.savefig('%s.png' % num_sites_graph_filename)
        P.close()

    @pylab_ioff
    def _plot_coverage_vs_E_value(self, start_results):
        """
        Plot how coverage of the start instances changes for the final discretised E-value
        """
        P.figure()
        P.scatter(
            [100 * start_result.em_result.coverage for start_result in start_results],
            [start_result.log_E_value for start_result in start_results],
        )
        P.xlim((-2, 102))
        P.xlabel('\% coverage of start instances')
        P.ylabel('log E-value')
        filename = os.path.join(
            self._output_dir_for_motif(), 'coverage-vs-E-value')
        logger.info('Saving coverage vs. E-value graph to %s', filename)
        P.savefig('%s.eps' % filename)
        P.savefig('%s.png' % filename)
        P.close()

    def __call__(self, starting_point_finder):
        """
        Find the instances by running EM on each start.

        Takes a list of instances already found that is updated.
        """
        #
        # for each start, run EM
        #
        level = logging.DEBUG
        start_results = []
        for starts_for_width in starting_point_finder.starting_points:
            for starts in starts_for_width:
                for start in starts.value:
                    # run EM algorithm
                    StartingPointFinder.log_start(
                        start, self.options.max_w, tag='starting EM', level=logging.DEBUG)
                    em_result = self._run_em_from_start(start)
                    sig = self.significance.log_E_value(em_result.model)
                    last_sig = None
                    em_results = starting_point_finder.em_results[
                        (start.seed, start.num_sites)]
                    if em_results:
                        last_sig = em_results[-1]
                        logger.log(
                            level, 'Significance has changed from %f to %f since last evaluation', last_sig, sig)
                    em_results.append(sig)
                    start_results.append(StartResult(sig, start, em_result))
                    #
                    # We don't want to run EM on all the starts we found as that can take too long.
                    # Here we decide whether to carry on with the rest of the starts for this width and number of sites.
                    # if we never ran EM on this start before or if the E-value stayed the same or improved
                    # since the last run of EM then we don't examine any more starts
                    #
                    if None == last_sig or sig <= last_sig:
                        break

        logger.info('Ran EM on %d starts', len(start_results))
        if not start_results:
            logger.warning(
                "Could not find any starts in the data for these parameter settings. "
                "You could try increasing the '--starts-per-motif' option"
            )
            return None

        #
        # Plot some stats about the starts
        #
        self._plot_num_sites_graph(start_results)
        self._plot_coverage_vs_E_value(start_results)

        #
        # sort the motifs to find the best one and make a logo for it.
        #
        start_results.sort()
        best = start_results[0]
        best_model = best.em_result.model
        consensus = consensus_from_pssm(
            N.exp(best_model.bs.pssm.log_probs.values()), self.bg_freqs)
        num_sites = int(best_model.bs.pssm.num_samples)
        logo(
            N.exp(best_model.bs.pssm.log_probs.values()),
            'STEME-motif-%02d' % self.motif_idx,
            self.options.output_dir,
            make_png=True,
            write_title=False
        )

        #
        # Create log odds matrix for p-values
        #
        scaled_log_odds, pssm_cdf = create_log_odds(
            self.bg_freqs,
            N.exp(best_model.bs.pssm.log_probs.values()),
            roundto=3
        )

        #
        # print the sites used to make the motif if requested
        #
        if self.options.print_sites:
            motif_sites_filename = os.path.join(
                self._output_dir_for_motif(), 'motif-sites.txt')
            with open(motif_sites_filename, 'w') as motif_sites_file:
                _avg_distance, _predictions = self._print_instances(
                    motif_sites_file, best.em_result.discretise_instances, best_model, scaled_log_odds, pssm_cdf, consensus)

        #
        # make predictions from best motif
        #
        instances = self._make_predictions_using_instance_finder(best_model)
        input_stats = analyse_num_instances(instances, self.input_sequences)
        logger.info(
            'Making predictions: log10 E-value=%.2f; consensus=%s, # sites used to generate motif=%d, predicted sites=%d',
            best.log_E_value /
            N.log(10.), consensus, num_sites, input_stats.count
        )
        logger.info(
            'Found %d instances in %d bases in %d INPUT      sequences for %.1e instances/base',
            input_stats.count, input_stats.num_bases, input_stats.num_seqs, input_stats.per_base
        )

        #
        # If given background data, also make predictions on that
        #
        bg_stats = None
        if self.bg_sequences:
            bg_prediction_model = self.bg_sequences.create_model(
                self.bg_bg_manager.get_bg_model(best_model.W), best_model.W)
            bg_prediction_model.bs.pssm.log_probs.values()[
                :] = best_model.bs.pssm.log_probs.values()
            bg_prediction_model.bs.recalculate()
            bg_prediction_model.lambda_ = best_model.lambda_
            instance_finder = FindInstances(
                self.bg_sequences.data, bg_prediction_model, self.options.prediction_Z_threshold)
            instance_finder()
            bg_stats = analyse_num_instances(
                instance_finder.instances, self.bg_sequences)
            logger.info(
                'Found %d instances in %d bases in %d BACKGROUND sequences for %.1e instances/base%s',
                bg_stats.count, bg_stats.num_bases, bg_stats.num_seqs, bg_stats.per_base,
                bg_stats.per_base and (
                    ' for a %.1f-fold enrichment in the input' %
                    (input_stats.per_base / bg_stats.per_base)) or ''
            )

        #
        # For each W-mer we found
        #
        predicted_sites_filename = os.path.join(
            self._output_dir_for_motif(), 'predicted-sites.txt')
        with open(predicted_sites_filename, 'w') as predicted_sites_file:
            avg_distance, predictions = self._print_instances(
                predicted_sites_file, instances, best_model, scaled_log_odds, pssm_cdf, consensus)

        #
        # Update the width-specific priors
        #
        self.input_sequences.data.update_priors()

        #
        # Print stats over all W-mers
        #
        logger.info('Average distance per base: %.3f', avg_distance)

        #
        # remove those starts that overlap the best start's original instances or
        # the instances that EM found
        #
        starting_point_finder.remove_overlapping(
            best.start.best_w_mers, best.start.model.W)
        non_overlapping_instances = InstanceVec()
        non_overlapping_instances.extend(instances)
        erase_overlapping_from_score_sorted(
            non_overlapping_instances, best_model.W)
        sort_instances_by_position(non_overlapping_instances)
        starting_point_finder.remove_overlapping(
            non_overlapping_instances, best.start.model.W)

        #
        # plot log likelihoods
        #
        with pylab_context_ioff():
            P.figure()
            P.plot(best.em_result.LLs)
            P.ylabel('log likelihood')
            P.xlabel('iteration')
            P.savefig(os.path.join(self._output_dir_for_motif(), 'LLs.png'))
            P.close()

        # Don't retain all the start results unless explicitly asked
        # as they can occupy a lot of memory
        if not self.options.store_start_results:
            start_results = None

        return Motif(
            idx=self.motif_idx,
            model=best_model,
            num_sites=num_sites,
            LLR=-1.,
            log_E_value=best.log_E_value,
            start_results=start_results,
            predictions=predictions,
            input_stats=input_stats,
            bg_stats=bg_stats
        )

    def _print_instances(self, f, instances, model, scaled_log_odds, pssm_cdf, consensus):
        """
        Print the instances to the file.
        """
        predictions = []
        avg_distance = 0
        f.write(
            "# W-mer, p(binding), p-value, strand, Hamming distance from consensus, sequence, offset\n")
        for w_mer in instances:

            #
            # Get the W-mer from the position and orientation
            #
            W_mer = self.input_sequences.data.get_W_mer(
                model.W, w_mer.global_pos)
            if w_mer.rev_comp:
                W_mer = reverse_complement(W_mer)
            seq, offset = self.input_sequences.data.pos_localise(
                w_mer.global_pos)
            interval = Interval(offset, offset + model.W)

            #
            # Score the W-mer to make and store a prediction
            #
            score = score_pssm_on_seq(
                scaled_log_odds, imap(index_for_base, W_mer))
            p_value = pssm_cdf[score]
            prediction = Prediction(
                seq=self.ids[seq],
                interval=interval,
                rev_comp=w_mer.rev_comp,
                p_value=p_value
            )
            predictions.append(prediction)

            #
            # Calculate and print some stats
            #
            distance = hamming_distance(W_mer, consensus)
            avg_distance += distance
            logger.debug(
                'Predicted site: %s score=%.5f p-value=%.1e strand=%s distance=%d seq=%d offset=%d',
                W_mer, w_mer.Z, p_value, w_mer.rev_comp and '-' or '+', distance, seq, offset
            )
            f.write(
                "%s,%.5f,%.1e,%s,%d,%d,%d\n" % (
                    W_mer, w_mer.Z, p_value, w_mer.rev_comp and '-' or '+', distance, seq, offset
                )
            )

            #
            # Adjust the position-specific prior to account for this binding site. I.e.
            # erase the binding site.
            #
            self.input_sequences.data.base_prior.adjust_for_binding_site(
                seq, offset, model.W, 1.)

        if len(instances):
            avg_distance /= float(len(instances) * model.W)

        return avg_distance, predictions


def _bg_model_creator_fn(options):
    """Returns a function to make a background Markov model for the order specified in the options.
    """
    try:
        logger.info('Using a Markov model of order %d' %
                    options.bg_model_order)
        return globals()['create_markov_model_order_from_index_%d' % options.bg_model_order]
    except AttributeError:
        raise ValueError(
            'STEME does not support Markov background models of order %d' %
            options.bg_model_order)


def _create_bg_model_from_Markov_fn(options):
    """Returns a function to make a background Markov model for the order specified in the options.
    """
    try:
        logger.info(
            'Creating a background model from a Markov model of order %d' %
            options.bg_model_order)
        return globals()['create_bg_model_from_Markov_model_%d' % options.bg_model_order]
    except AttributeError:
        raise ValueError(
            'STEME does not support Markov background models of order %d' %
            options.bg_model_order)


def build_Markov_model_of(index, options):
    """Build a Markov model of the sequences.
    """
    return _bg_model_creator_fn(options)(index, options.back_dist_prior)


class BgModelManager(object):

    """Manages background models.
    """

    def __init__(self, sequences, mm, options):
        self.sequences = sequences
        self.mm = mm
        self.options = options


class LogLikelihoodBgModelManager(BgModelManager):

    """Creates background models for various widths based on log likelihoods.
    """

    def __init__(self, sequences, mm, options):
        super(LogLikelihoodBgModelManager, self).__init__(
            sequences, mm, options)
        logger.info(
            'Calculating likelihoods of %d bases.', self.sequences.data.N)
        self.lls = self.mm.calculate_likelihoods(self.sequences.data)
        self.bg_models = {}

    def get_bg_model(self, W):
        """
        Return a background model of the sequences for the specified width.
        """
        if W not in self.bg_models:
            logger.info("Creating background model for W=%d", W)
            self.bg_models[W] = create_bg_model_from_base_likelihoods(
                W, self.sequences.data,
                self.lls,
                self.sequences.freqs_with_pseudo_counts
            )
        return self.bg_models[W]


class MarkovBgModelManager(BgModelManager):

    """Creates background models using a Markov model.
    """

    def __init__(self, sequences, mm, options):
        super(MarkovBgModelManager, self).__init__(sequences, mm, options)
        # we can use same background model for all widths
        self.bg_model = _create_bg_model_from_Markov_fn(options)(self.mm)

    def get_bg_model(self, W):
        """
        Return a background model of the sequences for the specified width.
        """
        return self.bg_model


class SequenceSet(object):

    """
    A set of sequences and associated useful data such as background models and an index.
    """

    def __init__(self, fasta_filename, options):
        self.options = options
        self.num_bases, self.seqs, self.ids, self.index = read_sequences(
            fasta_filename, self.options)
        logging.info('max_W=%d', self.options.max_w)
        self.data = Data(self.index, max_W=self.options.max_w)
        self.occs = occurrences_from_index(self.data.index)
        self.freqs = ZeroOrderFrequencies(list(self.occs[:4]))
        self.freqs_with_pseudo_counts = self.freqs.add_pseudo_counts(
            self.options.back_dist_prior)
        self.mm, _ = build_Markov_model_of(self.data.index, self.options)

    def composition(self):
        """
        The composition of the sequences.
        """
        freqs = self.occs[:4].copy()
        freqs /= freqs.sum()
        unknown = self.occs[4] / self.occs.sum()
        return self.num_bases, len(self.seqs), freqs, unknown

    def create_model(self, bg, W):
        "@return: A model."
        bs = PssmBindingSiteModel(
            initialise_uniform_pssm(W, self.options.alphabet_size))
        model = Model(self.data, bs, bg, _lambda=self.options.lambda_)
        return model

    def index_for_seq_id(self, seq_id):
        ids = [i for i, _id in enumerate(self.ids) if _id == seq_id]
        if not ids:
            raise ValueError('Could not find sequence with id = %s' % seq_id)
        if len(ids) > 1:
            raise ValueError(
                'Found %d sequences with id = %s' % (len(ids), seq_id))
        return ids[0]

    def get_bg_freqs(self, with_pseudo_counts=True, asarray=True):
        result = self.occs[:4].copy()
        if with_pseudo_counts:
            result += self.options.back_dist_prior
        result /= result.sum()
        if asarray:
            return result
        else:
            return list(result)

    def write_0_order_bg_freqs(self, filename):
        """
        Write the background frequencies as a 0-order Markov model in MEME format.
        """
        bg_freqs = self.get_bg_freqs()
        f = open(filename, 'w')
        print >>f, """# order 0
A %.3e
C %.3e
G %.3e
T %.3e""" % (bg_freqs[0], bg_freqs[1], bg_freqs[2], bg_freqs[3])


class Algorithm(object):

    """
    An object that applies the STEME algorithm to the sequences in a FASTA file.
    """

    @staticmethod
    def version():
        "@return: A tuple containing 3 version numbers."
        return __major_version__, __minor_version__, __release_version__

    def __init__(self, options):
        "Initialise."
        self.options = options
        if self.options.width:
            # make sure max_w is correct
            self.options.max_w = max(self.options.width)
        ensure_dir_exists(self.options.output_dir)
        self.name = 'STEM'
        self._setup_logging()
        self.bg_models = dict()
        self.outputs = []

    def create_model_of_input(self, W):
        return self.input_sequences.create_model(self.bg_manager.get_bg_model(W), W)

    def _setup_logging(self, level=logging.INFO):
        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s")
        file_handler = logging.FileHandler(
            os.path.join(self.options.output_dir, 'STEME.log'), mode='w')
        file_handler.setFormatter(formatter)
        file_handler.setLevel(level)
        logging.getLogger('').addHandler(file_handler)
        # uncomment following 3 lines to detect where logging format errors are raised.
#        def handleError(record):
#            raise RuntimeError(record)
#        file_handler.handleError = handleError

    def _set_min_and_max_number_sites(self):
        "Configure our max and min number of sites options if user has not set them explicitly."
        if not self.options.max_num_sites:
            # By default we will look for at most 3 * # seqs / 2 sites with a
            # minimum of 2
            self.options.max_num_sites = max(
                3 * len(self.input_sequences.seqs) / 2, 2)

            # Check if we've exceeded the hard limit
            hard_limit = 10000
            if self.options.max_num_sites > hard_limit:
                logger.warning(
                    'Implicit maximum number of sites (%d) is larger than hard limit (%d), capping at hard limit.', self.options.max_num_sites, hard_limit)
                logger.warning(
                    'Hit cap of %d on maximum number of sites, user can over-ride this with the --max-sites option. High values will slow STEME down.' % hard_limit)
                self.options.max_num_sites = hard_limit

            # however if user has set minimum number of sites always have at
            # least double this
            if self.options.min_num_sites:
                twice = 2 * self.options.min_num_sites
                if self.options.max_num_sites < twice:
                    self.options.max_num_sites = twice
                    logger.info(
                        'Maximum number of sites not explicitly set and is less than twice explicit minimum. Increasing to %d', self.options.max_num_sites)

            logger.info(
                'Maximum number of sites not explicitly set. Was set to %d.', self.options.max_num_sites)

        if not self.options.min_num_sites:
            # By default we will look for at least one site in every 10th
            # sequence
            self.options.min_num_sites = len(self.input_sequences.seqs) / 10

            # but always go down to at least half the maximum number of sites
            if self.options.min_num_sites > self.options.max_num_sites / 2:
                logger.info(
                    'Reseting minimum number of sites to half of maximum.')
                self.options.min_num_sites = self.options.max_num_sites / 2

            # looking for less than 2 sites makes no sense
            if self.options.min_num_sites < 2:
                logger.info(
                    'Makes no sense to have minimum number of sites less than two.')
                self.options.min_num_sites = 2

            logger.info(
                'Minimum number of sites not explicitly set. Was set to %d.', self.options.min_num_sites)

        if self.options.max_num_sites < 2:
            raise RuntimeError(
                'Bad configuration option max-num-sites: Makes no sense to look for a motif formed of less than 2 sites.')

        if self.options.min_num_sites < 2:
            raise RuntimeError(
                'Bad configuration option min-num-sites: Makes no sense to look for a motif formed of less than 2 sites.')

        if self.options.min_num_sites > self.options.max_num_sites:
            raise RuntimeError(
                'Bad configuration options: Minimum number of sites should be less than maximum number of sites.')

    def _initialise_bg(self):
        "Initialise a Markov model for the background."

        self.bg_sequences = None
        # Was a background FASTA file specified?
        if self.options.bg_fasta_file:
            # yes
            logger.info(
                'Building background Markov model from sequences in %s', self.options.bg_fasta_file)
            self.bg_sequences = SequenceSet(
                self.options.bg_fasta_file, self.options)
            self.bg_manager = get_background_manager(
                self.input_sequences, self.bg_sequences.mm, self.options)
            self.bg_bg_manager = get_background_manager(
                self.bg_sequences, self.bg_sequences.mm, self.options)
            self.freqs = self.bg_sequences.freqs
            self.freqs_with_pseudo_counts = self.bg_sequences.freqs_with_pseudo_counts
        else:
            # no, build background model from sequences we are searching for
            # motifs in
            logger.info(
                'Building background Markov model from input sequences')
            self.bg_manager = get_background_manager(
                self.input_sequences, self.input_sequences.mm, self.options)
            self.bg_bg_manager = None
            self.freqs = self.input_sequences.freqs
            self.freqs_with_pseudo_counts = self.input_sequences.freqs_with_pseudo_counts

    def background_composition(self):
        """
        The composition of the background sequences.
        """
        if self.bg_sequences:
            return self.bg_sequences.composition()
        else:
            return None, None, None, None

    def _initialise_p_value_tables(self):
        "Initialise MEME log likelihood ratio p-value tables."
        bg_freqs = [self.freqs_with_pseudo_counts.freq(b) for b in xrange(4)]
        # self.p_value_calculator = create_bejerano_pvalue_calculator(bg_freqs, self.options.max_num_sites)
        Q = self.options.pvalue_table_lattice_size and self.options.pvalue_table_lattice_size or 2 * \
            self.options.max_num_sites
        self.p_value_calculator = create_shifted_hirji_pvalue_calculator(
            bg_freqs, self.options.max_num_sites, Q)
        self.significance = Significance(
            self.input_sequences.data,
            self.input_sequences.freqs_with_pseudo_counts,
            self.p_value_calculator
        )

    def create_motif_finder(self, motif_idx=0):
        "@return: A MotifFinder object"
        return MotifFinder(
            self.input_sequences,
            self.bg_manager,
            self.significance,
            self.input_sequences.ids,
            motif_idx,
            self.input_sequences.get_bg_freqs(),
            self.options,
            self.bg_sequences,
            self.bg_bg_manager
        )

    def _check_options(self):
        if self.options.min_num_sites and self.options.max_num_sites and self.options.min_num_sites > self.options.max_num_sites:
            raise ValueError('Minimum number of sites > maximum : %d > %d' %
                             (self.options.min_num_sites,  self.options.max_num_sites))

    def _initialise(self, fasta):
        "Load sequences, calculate background frequencies and initialise p-value tables."
        self.timings = []
        with Timer(msg='initialise STEME') as timer:
            #
            # Check options
            #
            self._check_options()

            #
            # Load the sequences
            #
            with Timer(msg='load sequences') as _timer:
                self.fasta = fasta
                self.input_sequences = SequenceSet(fasta, self.options)
                self._set_min_and_max_number_sites()

            #
            # Calculate the background frequencies and likelihoods
            #
            with Timer(msg='initialise background model') as _timer:
                self._initialise_bg()

            #
            # Initialise the p-value tables
            #
            with Timer(
                msg='initialise p-value tables for numbers of sites from %d to %d' % (
                    self.options.min_num_sites, self.options.max_num_sites
                )
            ) as _timer:
                self._initialise_p_value_tables()

            #
            # Initialise outputs
            #
            self._initialise_minimal_meme_output()
            self._initialise_meme_like_output()
            self._initialise_html_output()

            #
            # Save timings
            #
            self.timings.append((timer._msg, timer.timer.duration()))

    def _initialise_meme_like_output(self):
        "Initialise MEME-like output if required"
        if self.options.meme_like_output:
            self.meme_like_output_file = os.path.join(
                self.options.output_dir, self.options.meme_like_output)
            output = MemeLikeOutput(open(self.meme_like_output_file, 'w'))
            output.initialise(self)
            self.outputs.append(output)

    def _initialise_minimal_meme_output(self):
        "Initialise minimal MEME output"
        self.minimal_meme_output_file = os.path.join(
            self.options.output_dir, 'steme.txt')
        output = MinimalMemeOutput(open(self.minimal_meme_output_file, 'w'))
        output.initialise(self)
        self.outputs.append(output)

    def _initialise_html_output(self):
        "Initialise HTML output if required and we can import it."
        if self.options.html_output:
            filename = os.path.join(
                self.options.output_dir, self.options.html_output)
            try:
                from ._html_output import HTMLOutput
            except ImportError as e:
                logger.warning('Could not import HTML output module: %s', e)
                logger.warning('Not writing HTML output to %s', filename)
            else:
                output = HTMLOutput(self.options.output_dir)
                output.initialise(self)
                self.html_output_file = filename
                self.outputs.append(output)

    def output_motif(self, motif, seconds_taken):
        "Write a motif to outputs."
        for output in self.outputs:
            output.found_motif(self, motif, seconds_taken=seconds_taken)

    def _finalise(self):
        "Terminate outputs and do post-processing."
        #
        # Finalise outputs
        #
        for output in self.outputs:
            output.finalise(self)

        #
        # Call TomTom if asked to
        #
        self._tomtom()

        #
        # Report how long various parts took
        #
        for msg, duration in self.timings:
            logging.info('%20s took %9.1f seconds (%2.0f%% of total)',
                         msg, duration, 100 * duration / self.total_duration)

    def _find_starts(self):
        "Find the starts we will use."
        with Timer(msg='find starts') as timer:
            self.start_finder = StartingPointFinder(
                self.input_sequences, self.bg_manager, self.significance, self.options)
            self.start_finder.find_starts()
            self.timings.append((timer._msg, timer.timer.duration()))

    def _find_motifs(self):
        "Find the motifs."
        with Timer(msg='find motifs') as timer:
            self.motifs = []
            for _motif_idx in xrange(self.options.num_motifs):
                if not self._find_motif():
                    # We couldn't find a motif, perhaps we ran out of starts
                    break
            self.timings.append((timer._msg, timer.timer.duration()))

    def _find_motif(self):
        "Find a motif."
        motif_idx = len(self.motifs)
        with Timer(msg='find motif %d' % motif_idx) as timer:
            #
            # Create and store Motif object
            #
            motif_finder = self.create_motif_finder(motif_idx)
            motif = motif_finder(self.start_finder)
            if None != motif:
                self.motifs.append(motif)
                self.output_motif(motif, seconds_taken=timer.timer.duration())
            return motif

    def _tomtom(self):
        """
        Run TOMTOM on the MEME output.
        """
        if self.options.tomtom and self.options.meme_like_output:
            logger.info(
                'Running TOMTOM tool on %s', ' '.join(self.options.tomtom))
            self.tomtom_dir = os.path.join(self.options.output_dir, 'TOMTOM')
            ensure_dir_exists(self.tomtom_dir)
            bg_filename = os.path.join(
                self.options.output_dir, 'bg-0-order.model')
            self.input_sequences.write_0_order_bg_freqs(bg_filename)
            args = [
                "tomtom",
                "-bfile", bg_filename,
                "-oc", self.tomtom_dir,
                self.meme_like_output_file
            ]
            # add the databases to the arg list
            args.extend(self.options.tomtom)
            logger.info('Calling TOMTOM with args: %s', ' '.join(args))
            tomtom_out = open(
                os.path.join(self.tomtom_dir, 'tomtom.stdout'), 'w')
            tomtom_err = open(
                os.path.join(self.tomtom_dir, 'tomtom.stderr'), 'w')
            subprocess.check_call(args, stdout=tomtom_out, stderr=tomtom_err)

    def __call__(self, fasta):
        "Run the method."
        logger.info('Running STEME %s', version_string())
        with Timer(msg='run STEME') as _timer:
            self._initialise(fasta)
            self._find_starts()
            self._find_motifs()
            self.total_duration = _timer.timer.duration()
            self._finalise()
