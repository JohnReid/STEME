#
# Copyright John Reid 2009, 2010, 2012
#

"""
Code to run MEME algorithm.
"""


import subprocess
import re
import time
import warnings
from cookbook.interval import Interval
from stempy import ensure_dir_exists, logging, os, parse_options
from cookbook.named_tuple import namedtuple


logger = logging.getLogger(__name__)

Start = namedtuple(
    'Start', 'w0 nsites0 cons0 cons w nsites sig em_time niters cons_after_em')


if False:
    warnings.warn('Using debug MEME')
    _meme_binary = '/home/john/local/debug/bin/meme.bin'
else:
    _meme_binary = '/home/john/local/bin/meme.bin'


def run_meme(fasta, options, extra_args=None):
    """
    Runs MEME.
    """

    # set up command line
    meme_cmd_args = [
        _meme_binary,
        fasta,
        '-maxsize', '100000000',
        '-oc', options.output_dir,
        '-dna',
        '-mod', 'anr',
        '-revcomp',
        '-nomatrim',
        #'-print_starts',
        '-nostatus',
        #'-print_all',
        #'-trace',
        #'-spmap', 'uni',
        #'-spfuzz', '0.5',
    ]
    if hasattr(options, 'max_iters') and options.max_iters:
        meme_cmd_args += ['-maxiter', str(options.max_iters)]
    if hasattr(options, 'W') and options.W:
        meme_cmd_args += ['-w', str(options.W)]
    if hasattr(options, 'min_w') and options.min_w:
        meme_cmd_args += ['-minw', str(options.min_w)]
    if hasattr(options, 'max_w') and options.max_w:
        meme_cmd_args += ['-maxw', str(options.max_w)]
    if hasattr(options, 'min_num_sites') and options.min_num_sites:
        meme_cmd_args += ['-minsites', str(options.min_num_sites)]
    if hasattr(options, 'max_num_sites') and options.max_num_sites:
        meme_cmd_args += ['-maxsites', str(options.max_num_sites)]
    if hasattr(options, 'num_motifs') and options.num_motifs:
        meme_cmd_args += ['-nmotifs', str(options.num_motifs)]
    if extra_args:
        meme_cmd_args += extra_args

    # rune MEME
    logger.info('Running MEME: %s', ' '.join(meme_cmd_args))
    popen = subprocess.Popen(
        meme_cmd_args,
        stdin=None, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=False)
    stdoutdata, stderrdata = popen.communicate()

    # save output
    stdout_filename = os.path.join(options.output_dir, 'meme.out')
    stderr_filename = os.path.join(options.output_dir, 'meme.err')
    print >> open(stdout_filename, 'w'), stdoutdata
    print >> open(stderr_filename, 'w'), stderrdata

    # check return code
    if popen.returncode:
        logger.warning(stdoutdata)
        logger.error(stderrdata)
        raise RuntimeError('Could not run MEME: retcode=%d' % popen.returncode)

    # look for start points in stdout
    #component_re = re.compile('component\s+(\d+): lambda=\s+([0-9.]+) ps=\s+([0-9.]+) cons0=(\w+) score=([0-9.]+)')
    start_re = re.compile(
        '\(start\)\s+([0-9]+)\s+([0-9.]+)\s+(\w+)\s+-->\s+(\w+)\s+w\s+([0-9]+)\s+nsites\s+([0-9.]+)\s+sig\s+([0-9.]+|inf)')
    em_time_re = re.compile(
        'EM \(without discretization\) took ([0-9.]+) seconds for ([0-9]+) iterations. Final consensus: ([ACGTX]+)')
    theta_re = re.compile(
        '\s+(0\.[0-9]+)\s+(0\.[0-9]+)\s+(0\.[0-9]+)\s+(0\.[0-9]+)')
    lambda_re = re.compile('lambda= (0\.[0-9]+) \(theta and obs\)')
    starts = []
    Zs = []
    thetas = []
    theta = None
    lambdas = []
    Z = None
    for l in stdoutdata.splitlines():
        if l.startswith('EM (without discretization) took'):
            match = em_time_re.match(l)
            if not match:
                raise ValueError('Could not parse EM timing line: %s', l)
            logger.debug('Got EM timing line from MEME output file\n%s', l)
            em_time = float(match.group(1))
            niters = int(match.group(2))
            cons_after_em = match.group(3)
            logging.info(cons_after_em)

#        if l.startswith('component'):
#            match = component_re.match(l)
#            if not match:
#                raise ValueError('Could not parse start line: %s', l)
#            logger.debug('Got MEME start point from MEME output file\n%s', l)
# w0, num_sites, start_cons, model_cons, w, num_sites_dis, sig = match.groups()
# model_cons, start_cons = match.groups()
#            component, lambda_, ps, start_cons, start_score = match.groups()
#            component = int(component)
#            ps = float(ps)
#            lambda_ = float(lambda_)
#            start_score = float(start_score)
# w0 = int(w0)
# num_sites = float(num_sites)
# w = int(w)
# num_sites_dis = int(num_sites_dis)
#            start = Start(component=component, seed=start_cons, ps=ps, lambda_=lambda_, score=start_score)
#            logger.info('MEME start: %s; num_sites=%.1f', start, start.ps * start.lambda_)
#            starts.append(start)
        elif l.startswith('(start)'):
            match = start_re.match(l)
            if not match:
                raise ValueError('Could not parse start line: %s' % l)
            logger.debug('Got MEME EM info from MEME output file\n%s', l)
            w0 = int(match.group(1))
            nsites0 = float(match.group(2))
            cons0 = match.group(3).strip()
            cons = match.group(4).strip()
            w = int(match.group(5))
            nsites = float(match.group(6))
            sig = float(match.group(7))
            start = Start(
                w0=w0, nsites0=nsites0, cons0=cons0, cons=cons, w=w, nsites=nsites,
                sig=sig, em_time=em_time, niters=niters, cons_after_em=cons_after_em)
            logger.info('EM start: %s', start)
            starts.append(start)
            # print w0, nsites0, cons0, cons, w, nsites, sig
        elif l.startswith('Z_ij: '):
            if l.startswith('Z_ij: lambda'):
                Zs.append(dict())
            else:
                _, i, j, strand, z = l.split()
                i = int(i)
                j = int(j)
                z = float(z)
                Zs[-1][(i, j, strand)] = z
        elif l.startswith('obs:'):
            theta = list()
            thetas.append(theta)
        elif l.startswith(' 0.'):
            match = theta_re.match(l)
            if match:
                theta.append(map(float, match.groups()))
        elif l.startswith('lambda= '):
            match = lambda_re.match(l)
            if match:
                lambdas.append(float(match.group(1)))

#    if len(em_infos) != len(starts):
#        raise RuntimeError('Should have found same number of starts as EM infos.')
    # return details
    logger.debug('Finished running MEME')
    if len(Zs) != len(thetas):
        raise ValueError('Parsed %d Zs and %d thetas', len(Zs), len(thetas))
    if len(Zs) != len(lambdas):
        raise ValueError('Parsed %d Zs and %d lambdas', len(Zs), len(lambdas))
    return meme_cmd_args, stdoutdata, starts, Zs, thetas, lambdas


class Algorithm(object):

    "An object that applies the MEME algorithm to the sequences in a FASTA file."

    def __init__(self, options):
        "Initialise."
        self.options = options
        self.name = 'MEME'

    def __call__(self, fasta):
        "Run the method."
        start_time = time.time()

        ensure_dir_exists(self.options.output_dir)

        predictions = []

        # run MEME
        self.meme_cmd_args, self.stdoutdata, self.starts, self.Zs, self.thetas, self.lambdas = run_meme(
            fasta, self.options)

        # parse output
        from Bio import Motif
        for motif in Motif.parse(open(os.path.join(self.options.output_dir, 'meme.txt')), "MEME"):
            for instance in motif.instances:
                # MEME parser seems to count from 1, not 0
                start = instance.start - 1
                prediction = instance.sequence_name, Interval(
                    start, start + motif.length), instance.strand == '-'
                predictions.append(prediction)

        logger.info('MEME took %.1f seconds', time.time() - start_time)

        return predictions


def add_options(option_parser):
    "Add options for MEME to the option parser."
    pass


if '__main__' == __name__:

    #
    # Set up the logging and options
    #
    import sys
    options, args = parse_options(add_options)

    # fasta = '/home/john/Data/NTNU-TF-search-dataset/datasets/model_real/M00724.fas'
    # fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../fasta/T00759trimRM-test-x2.fa'))
    # fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../fasta/T00759-tiny.fa'))
    fasta = '/home/john/Data/Tompa-data-set/Real/hm22r.fasta'
    # fasta = '/home/john/Data/GappedPssms/apr-2009/T99006trimRM.fa'
    # fasta = '/home/john/Data/GappedPssms/apr-2009/T99004trimRM.fa'
    options.output_dir = os.path.abspath(os.path.join('output', 'MEME'))
    # del options.max_num_sites
    # del options.min_num_sites

    method = Algorithm(options)
    predictions = method(fasta)
    from Bio.Motif.Parsers.MEME import MEMEParser
    parser = MEMEParser()
    meme_output = parser.parse(
        open(os.path.join(options.output_dir, 'meme.txt')))

    # Sequence, #start(inclusive), #end(inclusive), #nucleotides
    # M00441-1-embl, 1531, 1543, ggccacgtcaccg
    f = sys.stdout
    print >> f, '#Sequence, start(inclusive), end(inclusive), nucleotides, p-value'
    for motif in meme_output.motifs:
        print 'Motif: E-value: %f' % motif.evalue
        for instance in motif.instances:
            print >> f, "%10s, %5d, %5d, %s, %g" % (
                instance.sequence_name,
                instance.start,
                instance.start + motif.length - 1,
                str(instance),
                instance.pvalue,
            )
