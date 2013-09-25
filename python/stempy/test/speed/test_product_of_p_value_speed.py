#
# Copyright John Reid 2010, 2011
#


"""
Code to test how quick calculation of the product of p-values is.
"""

from setup_environment import init_test_env
init_test_env(__file__)

import stempy
import logging, os, time
logging.basicConfig(level=logging.DEBUG)


options = stempy.get_default_options()
fasta_file = os.path.join(fasta_dir(), 'T00759-small.fa')
num_bases, seqs, ids, index = stempy.read_sequences(fasta_file, options)
stem = stempy.Stem(index)

# get the background frequencies
occs = stempy.occurrences_from_index(index)
zero_order_freqs = stempy.ZeroOrderFrequencies(list(occs[:4]))
zero_order_freqs_with_pseudo_counts = zero_order_freqs.add_pseudo_counts(options.back_dist_prior)

W = 12
pssm_storage = stempy.initialise_random_pssm(W=W, alphabet_size=4, alpha=.1)
bs_model = stempy.PssmBindingSiteModel(pssm_storage)
bg_model = stempy.Markov0OrderBackgroundModel(stem.data, zero_order_freqs_with_pseudo_counts, W=W)
model = stempy.Model(bs_model, bg_model)

ns = (2, 5, 10, 50, 100, 500, 1000) #, 400, 600, 1000, 2000)
times = []
E_value_repeats = 100000
for n in ns:
    # initialise MEME log likelihood ratio p-value tables
    stempy.reset_MEME_llr_pv_tables() # can leak memory! (probably not much though)
    start = time.time()
    stempy.init_MEME_llr_pv_tables(options.min_num_sites, n, zero_order_freqs_with_pseudo_counts)
    duration = time.time() - start
    logging.info('Initialising LLR p-value tables up to %d sites took %.3f seconds.', n, duration)
    times.append(duration)

    stempy.__google_profiler_start('pop-speed.prof')
    model.bs.pssm.num_samples = n
    _ = model.log_E_value # do one first just to make sure everything is initialised
    start = time.time()
    for i in xrange(E_value_repeats):
        _ = model.log_E_value
    E_value_duration = time.time() - start
    logging.info('Calculating E-value %d times for %d sites took %.3f seconds.', E_value_repeats, n, E_value_duration)
    stempy.__google_profiler_stop()


import pylab as P
P.figure()
P.semilogy(ns, times)
P.xlabel('number of sites')
P.ylabel('seconds')
P.savefig('p-value-times.eps')
P.savefig('p-value-times.png')
P.savefig('p-value-times.pdf')
P.show()


