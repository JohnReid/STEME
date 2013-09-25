#
# Copyright John Reid 2012
#


"""
Code to analyse distances (spacings) between pairs of occurrences of TFs.
"""

import logging
logger = logging.getLogger(__name__)

from cookbook.named_tuple import namedtuple
from scipy.special import gammaln, digamma
import pyicl
import numpy
import pylab
import os
from collections import defaultdict
from .scan import footprint

PairOccurrence = namedtuple('PairOccurrence', 'spacing seq pos strand')
Spacing = namedtuple(
    'Spacing', 'primary secondary same_strand upstream distance')


def add_max_distance_option(parser):
    parser.add_option(
        "-d",
        "--max-distance",
        type=int,
        default=30,
        help="Only look for occurrences of motifs up to MAX_DISTANCE base pairs apart",
        metavar="MAX_DISTANCE"
    )


def add_options(parser):
    """Add options for spacing functionality to the parser.
    """
    from . import scan
    scan.add_options(parser)
    add_max_distance_option(parser)


def spacing_idx(max_distance, distance, upstream, same_strand):
    """Calculates an index into a spacing array.
    
    Args:
        max_distance: The maximum distance represented in the spacing array.
        distance: The distance to calculate the index for (should be in [0, max_distance].
        upstream: True if the secondary occurrence upstream of the primary occurrence.
        same_strand: True if the secondary occurrence on the same strand as the primary occurrence.
    
    Returns:
        An index in [0, 4*(max_distance+1)].
    """
    offset = max_distance + 1
    if distance not in pyicl.IntInterval(0, offset):
        raise ValueError('Bad distance')
    assert offset in pyicl.IntInterval(0, 2 * offset)
    idx = distance + int(upstream) * offset + int(same_strand) * 2 * offset
    assert idx in pyicl.IntInterval(0, 4 * offset)
    return idx


def spacing_from_idx(max_distance, idx):
    """Return the spacing for this index.
    """
    offset = max_distance + 1
    if idx not in pyicl.IntInterval(0, 4 * offset):
        raise ValueError('Bad index')
    same_strand = idx >= (2 * offset)
    idx %= 2 * offset
    upstream = idx >= offset
    idx %= offset
    assert idx in pyicl.IntInterval(0, offset)
    return idx, upstream, same_strand


def _test_idx_fns(max_distance):
    for distance, upstream, same_strand in (
        (0, False, False),
        (0, True, False),
        (0, False, True),
        (0, True, True),
        (1, False, False),
        (1, True, False),
        (1, False, True),
        (1, True, True),
        (max_distance - 1, False, False),
        (max_distance - 1, True, False),
        (max_distance - 1, False, True),
        (max_distance - 1, True, True),
        (max_distance, False, False),
        (max_distance, True, False),
        (max_distance, False, True),
        (max_distance, True, True),
    ):
        assert (distance, upstream, same_strand) == spacing_from_idx(
            max_distance, spacing_idx(max_distance, distance, upstream, same_strand))


# Uncomment to test index functions
# _test_idx_fns(9)


def can_be_primary(max_distance, primary_footprint, secondary_len, seq_length):
    """Is the primary motif too close to the boundaries of the sequence to be considered?
    """
    if max_distance + secondary_len > primary_footprint.first:
        return False
    if seq_length - max_distance - secondary_len <= primary_footprint.last:
        return False
    return True


def ln_factorial(n):
    return gammaln(n + 1)


def calc_ln_n_choose(obs):
    """Calculates the log of n choose the observed counts.
    """
    return ln_factorial(obs.sum()) - ln_factorial(obs).sum()


def calc_multinomial_ln_likelihood_uniform_dist(obs):
    """Calculates the log likelihood of the observed counts under a uniform multinomial distribution.
    """
    n = obs.sum()
    if 0 == n:
        return .0
    return calc_ln_n_choose(obs) - n * numpy.log(len(obs))


def calc_ln_gamma_factor(obs, alpha):
    """Calculates a factor involved in the likelihood of a multinomial distribution with a Dirichlet prior.
    """
    return (
        gammaln(alpha.sum())
        - gammaln(alpha).sum()
        + gammaln(alpha + obs).sum()
        - gammaln(alpha.sum() + obs.sum())
    )


def calc_multinomial_ln_likelihood_dirichlet_prior(obs, alpha):
    """Calculates the log likelihood of the observed counts under a Dirichlet prior.
    """
    return calc_ln_gamma_factor(obs, alpha) + calc_ln_n_choose(obs)


# def calc_llr_statistic(obs, alpha):
#    """Calculates the ratio of evidence in favour of a Dirichlet prior
#    defined by the alphas over a uniform distribution.
#    """
# return (obs * (digamma(alpha) + numpy.log(4 * len(alpha)))).sum() -
# obs.sum() * digamma(alpha.sum())


def calc_llr_statistic(obs, alpha):
    """Calculates the ratio of evidence in favour of a Dirichlet prior
    defined by the alphas over a uniform distribution. 
    """
    N = obs.sum()
    alpha0 = alpha.sum()
    D = len(obs)
    return gammaln(alpha0) - gammaln(alpha0 + N) + (gammaln(alpha + obs) + gammaln(alpha)).sum() + N * numpy.log(D)


def plot_spacings(max_distance, spacing_counts):
    """Plot the spacing counts in four subplots.
    """
    pylab.figure()
    subplot = 1
    ymax = spacing_counts.max()
    for same_strand in (True, False):
        for upstream in (True, False):
            ax = pylab.subplot(2, 2, subplot)
            if not same_strand:
                pylab.xlabel(upstream and 'upstream' or 'downstream')
            if not upstream:
                x = numpy.arange(0, max_distance + 1)
                ax.yaxis.set_label_position('right')
                ax.yaxis.tick_right()
            else:
                pylab.ylabel(
                    same_strand and 'Same strand' or 'Opposite orientation')
                x = numpy.arange(-max_distance, 1)
            y = [spacing_counts[spacing_idx(max_distance, int(abs(distance)), upstream, same_strand)]
                 for distance in x]
            pylab.bar(x, y, align='center', width=.5, color='grey')
            pylab.ylim(0, ymax)
            pylab.xlim(x[0] - .5, x[-1] + .5)
            ax.yaxis.set_major_locator(pylab.MaxNLocator(integer=True))
            ax.xaxis.set_major_locator(pylab.MaxNLocator(integer=True))
            subplot = subplot + 1


def yield_pairs(occurrences, seq_infos, options):
    """Yield pairs of occurrences under a given maximum distance.
    """
    for i1, occ1 in enumerate(occurrences):  # first occurrence
        footprint1 = footprint(occ1)
        seq_length = seq_infos[occ1.seq].length

        for occ2 in occurrences[i1+1]:

            # only interested in occurrences on same sequence
            if occ2.seq != occ1.seq:
                break

            footprint2 = footprint(occ2)
            if not footprint1.disjoint(footprint2):
                continue  # not interested in overlapping occurrences

            distance = footprint2.distance(footprint1)
            if distance > options.max_distance:
                break  # only interested in pairs up to some distance apart
            assert distance in pyicl.IntInterval(0, options.max_distance + 1)

            yield seq_length, occ1, footprint1, occ2, footprint2, distance


def make_primary_counting_handler(spacings):
    """Returns a pair handling function that counts spacings.
    """
    def handle_primary(max_distance, primary, secondary, distance, upstream, same_strand):
        """Update the spacings array.
        """
        idx = spacing_idx(max_distance, distance, upstream, same_strand)
        spacings[(primary.motif, secondary.motif)][idx] += 1

    return handle_primary


def make_pair_handler_from_primary_handler(primary_handler, ignore_close_to_end, options):
    """Makes a handler that considers a pair of instances in both configurations (primary, secondary)
    and (secondary, primary).
    
    The handler can choose to ignore those pairs that couldn't be counted in pair statistics
    as they are too close to the beginning or end of a sequence.
    """
    def handle_pair(seq_length, occ1, footprint1, occ2, footprint2, distance):
        # are they on the same strand?
        same_strand = bool(occ1.strand == occ2.strand)

        # handle occ1 as the primary occurrence
        # check whether occ1 is too close to start or end to be a primary motif
        if not ignore_close_to_end or can_be_primary(options.max_distance, footprint1, footprint2.size, seq_length):
            primary_handler(
                options.max_distance, occ1, occ2, distance, occ1.strand == '-', same_strand)

        # don't handle occ2 as the primary occurrence if it is the same motif
        # as the secondary occurrence
        if occ1.motif != occ2.motif:
            # handle occ2 as the primary occurrence
            # check whether occ2 is too close to start or end to be a primary motif
            if not ignore_close_to_end or can_be_primary(options.max_distance, footprint2, footprint1.size, seq_length):
                primary_handler(
                    options.max_distance, occ2, occ1, distance, occ2.strand == '+', same_strand)

    return handle_pair


def handle_pairs(occurrences, seq_infos, handler, options):
    """Handle pairs of occurrences under a given maximum distance.
    """
    for i1, occ1 in enumerate(occurrences):  # first occurrence
        footprint1 = footprint(occ1)
        seq_length = seq_infos[occ1.seq].length

        for occ2 in occurrences[i1+1:]:

            # only interested in occurrences on same sequence
            if occ2.seq != occ1.seq:
                break

            # don't use yet - check which scripts call this.
            # they all need property in options
            # ignore homodimers if requested
#            if options.ignore_homodimers and occ1.motif == occ2.motif:
#                continue

            footprint2 = footprint(occ2)
            if not footprint1.disjoint(footprint2):
                continue  # not interested in overlapping occurrences

            distance = footprint2.distance(footprint1)
            if distance > options.max_distance:
                break  # only interested in pairs up to some distance apart
            assert distance in pyicl.IntInterval(0, options.max_distance + 1)

            handler(seq_length, occ1, footprint1, occ2, footprint2, distance)


def count_all_pairs(occurrences, seq_infos, ignore_close_to_end, options):
    """Count all pairs of motifs in the occurrences.
    """
    spacings = defaultdict(lambda: numpy.zeros(
        4 * (options.max_distance + 1), dtype=numpy.uint))
    primary_handler = make_primary_counting_handler(spacings)
    pair_handler = make_pair_handler_from_primary_handler(
        primary_handler, ignore_close_to_end, options)
    handle_pairs(occurrences, seq_infos, pair_handler, options)
    return spacings


def parse_spacings(f):
    """Parse a list of spacings. Returns a map from motif pairs to sets of spacings
    """
    result = defaultdict(list)
    for line in f:
        if line.startswith('#'):
            continue  # ignore comments
        primary, secondary, same_strand, upstream, distance = line.strip(
        ).split()
        distance = int(distance)
        if 'S' != same_strand and 'C' != same_strand:
            raise ValueError('Expecting "S" or "C" for same strand field')
        same_strand = 'S' == same_strand
        if 'U' != upstream and 'D' != upstream:
            raise ValueError('Expecting "U" or "D" for upstream field')
        upstream = 'U' == upstream
        result[primary, secondary].append(
            Spacing(primary=primary, secondary=secondary,
                    same_strand=same_strand, upstream=upstream, distance=distance)
        )
    return result


spacing_header = '%30s %-25s %s %2s %4s' % (
    'Primary', 'Secondary', 'Strand', 'Up', 'Dist')


def spacing_str(spacing):
    return '%30s %-30s %s %2s %4d' % (
        spacing.primary, spacing.secondary, spacing.same_strand and "S" or "C", spacing.upstream and "U" or "D", spacing.distance
    )


def check_spacings(motif_pairs, motifs):
    for primary, secondary in motif_pairs:
        if primary not in motifs:
            raise ValueError('Motif "%s" is not in the occurrences' % primary)
        if secondary not in motifs:
            raise ValueError(
                'Motif "%s" is not in the occurrences' % secondary)


def log_spacings(spacings):
    """Log the spacings
    """
    logger.info('Using following spacings:')
    logger.info(spacing_header)
    logger.info('*' * len(spacing_header))
    for pair_spacings in spacings:
        for spacing in pair_spacings:
            logger.info(spacing_str(spacing))


def load_spacings(spacings_filename, motifs, options):
    """Load spacings from filename. Check they are consistent with motif list.
    """
    logger.info('Reading spacings to find from: %s', spacings_filename)
    spacings = parse_spacings(open(spacings_filename))
    check_spacings(spacings.keys(), motifs)
    log_spacings(spacings.values())
    return spacings
