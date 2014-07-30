#
# Copyright John Reid 2012, 2013, 2014
#

"""
Code to analyse results of STEME PWM scanning.
"""


import logging
logger = logging.getLogger(__name__)

from . import html_copy_static
import os
import pyicl
import pylab
import numpy
import bisect
from cookbook.named_tuple import namedtuple
from collections import defaultdict
from itertools import ifilter
from cookbook.pylab_utils import pylab_context_ioff, create_format_cycler
# from cookbook.pylab_utils import simple_marker_styles


SeqInfo = namedtuple('SeqInfo', 'name length')
Occurrence = namedtuple(
    'Occurrence', 'motif wmer seq pos strand Z score pvalue')


def footprint(occ):
    """Return the footprint (interval) of the occurrence.
    """
    return pyicl.IntInterval(occ.pos, occ.pos + len(occ.wmer))


def parse_occurrence(line):
    """Parse one occurrence in the format outputted by steme-pwm-scan.
    """
    fields = line.strip().split(',')
    if 8 != len(fields):
        raise RuntimeError('Wrong number of fields in line: %s' % line)
    return Occurrence(
        motif=fields[0],
        wmer=fields[1],
        seq=int(fields[2]),
        pos=int(fields[3]),
        strand=fields[4],
        Z=float(fields[5]),
        score=float(fields[6]),
        pvalue=float(fields[7]),
    )


def line_is_not_comment(line):
    """Is the line a comment? Comments start with '#'
    """
    return not line.startswith('#')


def parse_occurrences(f):
    """Parse lines of occurrences in the format outputted by steme-pwm-scan.
    """
    return map(parse_occurrence, ifilter(line_is_not_comment, f))


def parse_seq_info(line):
    """Parse one sequence info in the format outputted by steme-pwm-scan.
    """
    fields = line.strip().split(',')
    if 2 != len(fields):
        raise RuntimeError('Wrong number of fields in line: %s' % line)
    return SeqInfo(
        name=fields[1],
        length=int(fields[0])
    )


def parse_seq_infos(f):
    """Parse lines of sequence infos in the format outputted by steme-pwm-scan.
    """
    return map(parse_seq_info, ifilter(line_is_not_comment, f))


def load_occurrences(options):
    """Load the occurrences and associated sequence lengths.
    """
    occurrences_filename = os.path.join(
        options.results_dir, 'steme-pwm-scan.out')
    seqs_filename = os.path.join(options.results_dir, 'steme-pwm-scan.seqs')
    logger.info(
        'Reading occurrences from: %s and sequence information from: %s',
        occurrences_filename, seqs_filename)
    return load_occurrences_from_stream(
        open(occurrences_filename), open(seqs_filename))


def load_occurrences_from_stream(occurrences_stream, seqs_stream):
    """Load the occurrences and associated sequence lengths.
    """
    #
    # Read in the occurrences
    #
    occurrences = parse_occurrences(occurrences_stream)

    #
    # Read in the sequence lengths
    #
    seq_infos = parse_seq_infos(seqs_stream)

    #
    # Gather all the motifs
    #
    motifs = set(occ.motif for occ in occurrences)

    #
    # Sort the occurrences by position
    #
    logger.info('Sorting %d occurrences', len(occurrences))
    occurrences.sort(key=lambda x: (x.seq, x.pos))

    return occurrences, seq_infos, motifs


def add_options(parser):
    """Add options for spacing functionality to the parser.
    """
    parser.add_option(
        "-r",
        "--results-dir",
        default='.',
        help="Look for results from PWM scan in DIR",
        metavar="DIR"
    )


def get_sites_by_motif(occurrences):
    """Organise occurrences by motif.
    """
    by_motif = defaultdict(list)
    for occ in occurrences:
        by_motif[occ.motif].append(occ)
    return by_motif


def plot_scores_per_motif(motifs, by_motif, format_cycler):
    """Plot the distribution of scores for each motif.
    """
    lines = []
    for i, motif in enumerate(motifs):
        occs = by_motif[motif]
        occs.sort(key=lambda occ: -occ.Z)
        fmt = format_cycler(i)
        lines.append(
            pylab.plot(
                numpy.arange(len(occs)),
                # numpy.linspace(0, 1, num=len(occs)),
                [occ.Z for occ in occs],
                label=motif,
                **fmt
            )[0]
        )
    pylab.ylim(ymax=1)
    pylab.ylabel('Z')
    pylab.xlabel('sites')
    pylab.title('Z by motif')
    return lines


def plot_scaled_positions(occs, seq_infos, **kwargs):
    scaled_positions = [(occ.pos + len(occ.wmer) / 2.)
                        / seq_infos[occ.seq].length for occ in occs]
    scaled_positions.sort()
    return pylab.plot(
        numpy.linspace(0, 1, num=len(scaled_positions)),
        scaled_positions,
        **kwargs
    )[0]


def plot_site_positions(motifs, occs, by_motif, seq_infos, format_cycler):
    """Plot the positions of the sites in the sequences.
    """
    lines = []
    lines.append(plot_scaled_positions(
        occs, seq_infos, label='ALL MOTIFS', c='k', linestyle='-'))
    for i, motif in enumerate(motifs):
        lines.append(plot_scaled_positions(
            by_motif[motif], seq_infos, label=motif, **format_cycler(i)))
    pylab.xlim(-.01, 1.01)
    pylab.gca().get_xaxis().set_visible(False)
    pylab.ylim(-.01, 1.01)
    pylab.gca().get_yaxis().set_ticks((0, .5, 1))
    pylab.gca().get_yaxis().set_ticklabels(('start', 'centre', 'end'))
    pylab.title('Positions in sequences')
    return lines


def calculate_num_sites_per_base(occs, seq_infos):
    """Calculate the density of sites in each sequence.
    """
    num_sites_per_seq = numpy.zeros(len(seq_infos), dtype=int)
    for occ in occs:
        num_sites_per_seq[occ.seq] += 1
    W = occs and len(occ.wmer) or 0
    return numpy.array([
        num_sites_per_seq[seq] and float(
            num_sites_per_seq[seq]) / (seq_info.length - W + 1) or 0
        for seq, seq_info
        in enumerate(seq_infos)
    ])


def plot_num_sites_per_seq(num_sites_per_base, **kwargs):
    return pylab.plot(
        numpy.arange(len(num_sites_per_base)),
        num_sites_per_base,
        **kwargs
    )


def adjust_sequence_xaxis(axes, seq_infos):
    axes.get_xaxis().set_ticks(())
    axes.get_xaxis().set_ticklabels(())
    pylab.xlabel('sequences')
    axis_adjust = len(seq_infos) / 100.
    pylab.xlim((- axis_adjust, len(seq_infos) - 1 + axis_adjust))


def calculate_per_motif_density(motifs, by_motif, seq_infos):
    """Calculate the sequence site density per motif.
    """
    num_sites_per_base = numpy.empty((len(motifs), len(seq_infos)))
    for m, motif in enumerate(motifs):
        num_sites_per_base[m] = calculate_num_sites_per_base(
            by_motif[motif], seq_infos)
    return num_sites_per_base


def calculate_motif_best_Z_per_sequence(motifs, by_motif, numseqs):
    """Calculate a numpy array indexed by motif, then sequence
    that represents the best score that motif had in that sequence
    """
    result = numpy.zeros((len(motifs), numseqs))
    for m, motif in enumerate(motifs):
        result_row = result[m]
        for occ in by_motif[motif]:
            assert occ.motif == motif
            result_row[occ.seq] = max(result_row[occ.seq], occ.Z)
    return result


def hier_cluster_and_permute(matrix):
    import scipy.cluster.hierarchy as hier
    from scipy.spatial.distance import pdist

    return hier.centroid(matrix)

    D = pdist(matrix)  # upper triangle of distance matrix as vector
    Y = hier.linkage(D, method='single')  # Cluster

    # return permuted matrix and dendrogram
    return Y


def plot_seq_coverage(best_Z, format_cycler):
    """Plot what proportion of sequences have sites for each motif at
    each Z-score threshold."""
    lines = []
    for i, motif_best_Z in enumerate(best_Z):
        sorted_best_Z = numpy.sort(motif_best_Z)
        first_non_zero = bisect.bisect(sorted_best_Z, 0)
        lines.append(pylab.plot(
            numpy.arange(best_Z.shape[1] - first_non_zero),
            sorted_best_Z[first_non_zero:][::-1],
            **format_cycler(i)))
    pylab.ylabel('Z')
    pylab.xlabel('sequences')
    return lines


def plot_seq_coverage_lines(best_Z, format_cycler):
    """Plot what proportion of sequences have sites for each motif at
    each Z-score threshold."""
    lines = []
    for i, motif_best_Z in enumerate(best_Z):
        sorted_best_Z = numpy.sort(motif_best_Z)
        first_non_zero = bisect.bisect(sorted_best_Z, 0)
        lines.append(pylab.plot(
            numpy.arange(best_Z.shape[1] - first_non_zero),
            sorted_best_Z[first_non_zero:][::-1],
            **format_cycler(i)))
    pylab.ylabel('Z')
    pylab.xlabel('sequences')
    return lines


def num_seq_clusters(num_seqs):
    """A heuristic to choose a number of clusters for the sequences based on
    how many motifs there are."""
    return max(2, int(numpy.log(num_seqs)))


def plot_best_Z(motifs, best_Z):
    """Plot the best Z for each motif in each sequence.
    """
    import scipy.cluster.hierarchy as hier
    import scipy.cluster.vq as vq
    fig = pylab.gcf()

    # Cluster (hiearchical) Y axis
    Y = hier.centroid(best_Z)
    axdendro = fig.add_axes([0.01, 0.02, 0.18, 0.96])
    axdendro.set_xticks([])
    axdendro.set_frame_on(False)
    dendro = hier.dendrogram(Y, labels=motifs, orientation='right')
    best_Z_permuted = best_Z[dendro['leaves'], :]

    # K-means cluster X axis
    xcentroid, xlabel = vq.kmeans2(
        best_Z.T, k=num_seq_clusters(best_Z.shape[1]))
    best_Z_permuted = best_Z_permuted[:, numpy.argsort(xlabel)]

    # Plot matrix
    axmatrix = fig.add_axes([0.4, 0.02, 0.5, 0.96])
    im = axmatrix.matshow(best_Z_permuted, aspect='auto', origin='lower')
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar
    axcolor = fig.add_axes([0.91, 0.02, 0.02, 0.96])
    pylab.colorbar(im, cax=axcolor)


def plot_collinearity(motifs, best_Z):
    """Plot the cooccurrences of motifs.
    """
    import scipy.cluster.hierarchy as hier
    # from scipy.stats import pearsonr
    M = len(motifs)
    cooccurrences = numpy.ones((M, M))
    for m1 in xrange(M):
        for m2 in xrange(M):
            # both = sum(numpy.logical_and(m1seqs, m2seqs))
            # cooccurrences[m1,m2] = both/float(sum(m2seqs))
            cooccurrences[m1, m2] = \
                numpy.sqrt(sum(best_Z[m1] * best_Z[m2])) \
                / numpy.linalg.norm(best_Z[m2])
            # rho, _ = pearsonr(best_Z[m1], best_Z[m2])
            # cooccurrences[m1, m2] = rho
    Y = hier.centroid(cooccurrences)
    index = hier.fcluster(Y, -1) - 1
    cooccurrences = cooccurrences[index, :]
    cooccurrences = cooccurrences[:, index]
    pylab.pcolor(cooccurrences)
    pylab.colorbar()
    ax = pylab.gca()
    ax.set_xticks([])
    # ax.set_xticks(.5 + numpy.arange(M))
    # ax.set_xticklabels(motifs)
    ax.set_yticks(.5 + numpy.arange(M))
    ax.set_yticklabels(numpy.asarray(motifs)[index])
    ax.set_xlim((0, M))
    ax.set_ylim((0, M))
    for line in ax.yaxis.get_ticklines():
        line.set_markersize(0)
    pylab.gcf().subplots_adjust(left=.27, bottom=.02, top=.98, right=.99)


def plot_seq_distribution(motifs, by_motif, seq_infos, format_cycler):
    """Plot the number of sites over the sequences.
    """
    pylab.title('Number of sites by sequence')
    lines = []
    num_sites_per_base = calculate_per_motif_density(
        motifs, by_motif, seq_infos)
    # overall_site_density = num_sites_per_base.sum(axis=0)
    # sortidx = overall_site_density.argsort()
    for m, motif in enumerate(motifs):
        num_sites_per_base[m].sort()
        lines.append(plot_num_sites_per_seq(
            num_sites_per_base[m, ::-1], label=motif, **format_cycler(m)))
    # pylab.gca().set_yscale('log')
    pylab.ylim(ymin=0)
    pylab.ylabel('sites per base')
    adjust_sequence_xaxis(pylab.gca(), seq_infos)
    return lines


def plot_seq_lengths(seq_infos):
    """Plot sequence lengths.
    """
    pylab.title('Sequence lengths')
    lengths = [info.length for info in seq_infos]
    pylab.scatter(
        numpy.arange(len(seq_infos)),
        lengths,
        alpha=.2,
        marker='o',
        c='k'
    )
    lengths.sort()
    pylab.plot(
        numpy.arange(len(seq_infos)),
        lengths,
        linestyle='-',
        linewidth=2,
        c='g'
    )
    adjust_sequence_xaxis(pylab.gca(), seq_infos)
    pylab.ylim(ymin=0)
    pylab.ylabel('base pairs')


def plot_occs_by_motif(by_motif):
    """Plot # occurrences for each motif.
    """
    sizes = [
        (len(occs), sum(occ.Z for occ in occs), name)
        for name, occs in by_motif.iteritems()]
    # expected = [(len(occs), name) for name, occs in by_motif.iteritems()]
    sizes.sort()
    bar_positions = numpy.arange(len(sizes))
    num_occs = numpy.asarray([s for s, e, n in sizes])
    total_Z = numpy.asarray([e for s, e, n in sizes])
    pylab.barh(
        bar_positions,
        num_occs,
        # left=total_Z,
        height=.8,
        align='center',
        label='Sites',
        color='blue',
    )
    pylab.barh(
        bar_positions,
        total_Z,
        height=.8,
        align='center',
        label='Total Z',
        color='blue',
        edgecolor='white',
        hatch='/',
    )
    pylab.yticks(bar_positions, [n for x, e, n in sizes])
    pylab.ylim(ymin=-.5, ymax=len(sizes) - .5)
    pylab.xlabel('occurrences')
    pylab.legend(loc='lower right')


def savefig(tag, options):
    """Save a figure to the results directory.
    """
    pylab.savefig(
        os.path.join(options.results_dir, 'scan-stats', '%s.png' % tag))


def create_figures(motifs, occs, by_motif, seq_infos, options):
    """Create figures.
    """

    from stempy import ensure_dir_exists
    ensure_dir_exists(os.path.join(options.results_dir, 'scan-stats'))

    # Size of figlegend
    if len(motifs) > 30:
        size = 6
    elif len(motifs) > 16:
        size = 8
    elif len(motifs) > 10:
        size = 10
    else:
        size = 12
    figlegendprops = {'size': size}

    # Format cycler for line plots
    format_cycler = create_format_cycler(
        linestyle=['--', '-.', '-', ':'],
        c=("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
           "#D55E00", "#CC79A7"))

    # Format cycler for marker plots
    # format_cycler_marker = create_format_cycler(
    #    marker=simple_marker_styles,
    #    c=("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    #       "#D55E00", "#CC79A7"))

    # Scan scores
    pylab.figure(figsize=(6, 4))
    lines = plot_scores_per_motif(motifs, by_motif, format_cycler)
    savefig('scan-scores', options)
    pylab.close()

    # Scan legend
    pylab.figure(figsize=(4.25, 4))
    pylab.figlegend(lines, motifs, 'center', prop=figlegendprops)
    savefig('scan-legend', options)
    pylab.close()

    # Best Z for each motif/sequence combination
    pylab.figure(figsize=(6, 4))
    best_Z = calculate_motif_best_Z_per_sequence(
        motifs, by_motif, len(seq_infos))
    plot_best_Z(motifs, best_Z)
    savefig('scan-best-Z', options)
    pylab.close()

    # Scan motif cooccurrences
    pylab.figure(figsize=(6, 4))
    # pylab.figlegend(lines, motifs, 'center')
    plot_collinearity(motifs, best_Z)
    savefig('scan-collinearity', options)
    pylab.close()

    # Scan positions
    pylab.figure(figsize=(6, 4))
    lines = plot_site_positions(motifs, occs, by_motif, seq_infos,
                                format_cycler)
    savefig('scan-positions', options)
    pylab.close()

    # Scan legend with all
    pylab.figure(figsize=(4.25, 4))
    pylab.figlegend(
        lines, ['ALL MOTIFS'] + motifs, 'center', prop=figlegendprops)
    savefig('scan-legend-with-all', options)
    pylab.close()

    # Sequence coverage
    pylab.figure(figsize=(6, 4))
    plot_seq_coverage(best_Z, format_cycler)
    savefig('scan-seq-coverage', options)
    pylab.close()

    # Scan sequences
    pylab.figure(figsize=(6, 4))
    lines = plot_seq_distribution(motifs, by_motif, seq_infos, format_cycler)
    savefig('scan-sequences', options)
    pylab.close()

    # Scan legend with markers
    # fig = pylab.figure(figsize=(4.25, 4))
    # pylab.figlegend(lines, motifs, 'center', prop=figlegendprops)
    # savefig('scan-legend-marker', options)
    # pylab.close()

    # Scan lengths
    pylab.figure(figsize=(6, 4))
    plot_seq_lengths(seq_infos)
    savefig('scan-lengths', options)
    pylab.close()

    # Scan occurrences by motif
    pylab.figure(figsize=(6, len(by_motif) / 4.))
    pylab.subplots_adjust(left=.3, bottom=.1, right=.96, top=.98)
    plot_occs_by_motif(by_motif)
    savefig('scan-occs-by-motif', options)
    pylab.close()


def create_html_output(dataset_name, motifs, occurrences, by_motif, seq_infos,
                       options):
    """Create HTML output.
    """
    from jinja2 import Environment, PackageLoader
    env = Environment(loader=PackageLoader('stempy', 'templates'))
    template = env.get_template('scan-stats.html')

    # copy the static info
    static_dir = os.path.join(options.results_dir, 'static')
    html_copy_static(static_dir)

    # write the HTML
    filename = os.path.join(options.results_dir, 'scan-stats.html')
    logger.info('Writing STEME scan statistics as HTML to %s', filename)
    num_bases = sum(info.length for info in seq_infos)
    with open(filename, 'w') as f:
        variables = {
            'dataset_name': dataset_name,
            'num_sites': len(occurrences),
            'num_motifs': len(motifs),
            'num_seqs': len(seq_infos),
            'num_bases': num_bases,
            'options': options,
            'num_seq_clusters': num_seq_clusters(len(seq_infos)),
        }
        f.write(template.render(**variables))

    # create the figures
    if len(occurrences):
        with pylab_context_ioff():
            create_figures(motifs, occurrences, by_motif, seq_infos, options)


def write_seq_centric_stats(out, motifs, occurrences, seq_infos, options):
    """Write sequence-centric stats in CSV format to a file."""
    def zero_num():
        return numpy.zeros(len(seq_infos), dtype=int)

    def zero_exp():
        return numpy.zeros(len(seq_infos), dtype=float)

    num_hits = defaultdict(zero_num)  # Number of hits
    exp_hits = defaultdict(zero_exp)  # Expected hits
    for occ in occurrences:
        num_hits[occ.motif][occ.seq] += 1
        exp_hits[occ.motif][occ.seq] += occ.Z
    motifs = num_hits.keys()
    print >>out, '# Length,ID,Total,Expected,%s' % ','.join(
        '%s,E(%s)' % (m, m) for m in motifs)
    for seq, seqinfo in enumerate(seq_infos):
        print >>out, "%d,%s,%d,%.4f,%s" % (
            seqinfo.length,
            seqinfo.name,
            sum(num_hits[m][seq] for m in motifs),
            sum(exp_hits[m][seq] for m in motifs),
            ','.join('%d,%.4f' % (num_hits[m][seq], exp_hits[m][seq])
                     for m in motifs))
