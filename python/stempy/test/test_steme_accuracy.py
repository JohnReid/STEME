#
# Copyright John Reid 2011, 2012
#

"""
Test STEME's performance on some small data sets of random sequences with planted sites.
"""

from setup_environment import init_test_env, logging, fasta_dir, is_debug_python
init_test_env(__file__, level=logging.INFO)

import stempy, os
from stempy.planted_sites import parse_fasta_for_sites, parse_meme_output_for_sites, calculate_positives_and_negatives
from infpy.roc import RocCalculator
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--run-meme", action='store_true', help="Run MEME as well on the sequences")
cmd_line_options, args = parser.parse_args()

#
# The data sets with minimum sensitivity and specificity values required.
#
if is_debug_python():
    logging.info('Detected debug version of Python, only using smallest data set.')
    data_sets = [
        #('random-seqs-03-050'        , .13, .86),
        ('random-seqs-with-Ns-03-050', .13, .86),
    ]
else:
    data_sets = [
        ('random-seqs-03-050'        , .13, .86),
        ('random-seqs-with-Ns-03-050', .13, .86),
        ('random-seqs-05-100'        , .40, .89),
        ('random-seqs-with-Ns-05-100', .60, .90),
        ('random-seqs-05-100'        , .40, .89), # cannot achieve (.6,.9) stats when finding starts up-front
        ('random-seqs-10-100'        , .60, .91),
        ('random-seqs-with-Ns-10-100', .20, .98), # lower specificity with Ns
        ('random-seqs-30-200'        , .46, .99),
        ('random-seqs-with-Ns-30-200', .46, .99),
    ]

rocs = dict()
meme_rocs = dict()
for data_set, min_sensitivity, min_specificity in data_sets:
    #
    # Set up the options
    #
    fasta_file = os.path.join(fasta_dir(), '%s.fasta' % data_set)
    options = stempy.get_default_options()
    options.min_w = 6
    options.max_w = 11
    options.output_dir = os.path.join('output', 'test-steme-accuracy', data_set)
    meme_output = os.path.join(options.output_dir, options.meme_like_output)
    
    #
    # Run the STEME algorithm
    #
    algorithm = stempy.Algorithm(options)
    algorithm(fasta_file)

    #
    # Analyse output
    #    
    records, true_sites = parse_fasta_for_sites(fasta_file)
    predicted_sites = parse_meme_output_for_sites(meme_output)
    TP, TN, FN, FP = calculate_positives_and_negatives(records, true_sites, predicted_sites)
    roc = RocCalculator(TP, FP, TN, FN)
    rocs[data_set] = roc
    
    #
    # Run MEME as well if asked to
    #
    if cmd_line_options.run_meme:
        from stempy import meme
        options.output_dir = os.path.join(options.output_dir, 'MEME')
        meme_algorithm = meme.Algorithm(options)
        meme_algorithm(fasta_file)
        
        #
        # Analyse MEME output
        #    
        records, true_sites = parse_fasta_for_sites(fasta_file)
        predicted_sites = parse_meme_output_for_sites(os.path.join(options.output_dir, 'meme.txt'))
        TP, TN, FN, FP = calculate_positives_and_negatives(records, true_sites, predicted_sites)
        meme_roc = RocCalculator(TP, FP, TN, FN)
        meme_rocs[data_set] = meme_roc
        logging.info("MEME's results: %s", meme_roc)
            
    #
    # Check we found enough sites and didn't predict too many
    #
    logging.info("STEME's results: %s", roc)
    if roc.specificity() < min_specificity:
        raise ValueError('Specificity too low: %s: %f < %f' % (data_set, roc.specificity(), min_specificity))
    if roc.sensitivity() < min_sensitivity:
        raise ValueError('Sensitivity too low: %s: %f < %f' % (data_set, roc.sensitivity(), min_sensitivity))

