#
# Copyright John Reid 2011
#


"""
Generate sequences according to simple background model. Also plant sites in sequences.
"""
from Bio import SeqIO
import logging
import numpy.random as R
import sys
from optparse import OptionParser
from stempy.planted_sites import seq_generator, plant_sites


logging.basicConfig(level=logging.INFO)

parser = OptionParser()
parser.add_option("-S", "--num-seqs", default=5,
                  type=int, help="Number of sequences")
parser.add_option("-N", "--avg-len", default=100,
                  type=int, help="Average length of sequences")
parser.add_option("-c", "--consensus", action="append",
                  help="Binding site consensus (can have more than one)")
parser.add_option("-p", "--p-base-changed", default=.2, type=float,
                  help="Probability any given binding site base differs from consensus")
parser.add_option("-o", "--output-file",
                  default="random-seqs.fasta", help="FASTA output file-name")
parser.add_option("-s", "--rng-seed", default=0,
                  type=int, help="Seed for random number generator")

sys.argv = [a.encode(sys.stdin.encoding or 'ascii') for a in sys.argv]
(options, args) = parser.parse_args()

#
# log options
#
for option in parser.option_list:
    if option.dest:
        logging.info('%32s: %-32s * %s', option.dest,
                     str(getattr(options, option.dest)), option.help)

#
# Seed RNG if asked to
#
if options.rng_seed:
    logging.info('Seeding random number generator with %d', options.rng_seed)
    R.seed(options.rng_seed)
else:
    logging.info('Not seeding random number generator with')


#
# Generate random sequences
#
logging.info('Generating %d sequences of average length %d',
             options.num_seqs, options.avg_len)
sequences = map(seq_generator(options.avg_len), xrange(options.num_seqs))


#
# Plant some sites
#
if not options.consensus:
    logging.warning('Was not given any consensus to plant in sequences.')
for consensus in options.consensus:
    plant_sites(sequences, consensus, options.p_base_changed)


#
# Write output
#
logging.info('Writing sequences to %s', options.output_file)
SeqIO.write(sequences, open(options.output_file, "w"), "fasta")
