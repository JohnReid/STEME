#
# Copyright John Reid 2011, 2012
#

"""
Test Markov background model.
"""

#
# Trickery to find update path to import stempy from
#
from setup_environment import init_test_env, logging
init_test_env(__file__)

import stempy
from math import exp, fabs
from itertools import chain, starmap

def pairs(iterable):
    """@return: Yield all pairs of consecutive elements of the iterable."""
    at_start = True
    last = None
    for v in iterable:
        if at_start:
            at_start = False
        else:
            yield last, v
        last = v

def base_ll(last, this):
    return this - last

def base_lls(lls):
    "@return: Yield log likelihood of each base."
    return starmap(base_ll, pairs(chain((0,), lls)))

def create_data(*seqs):
    string_set = stempy.StringSet()
    for seq in seqs:
        logging.info('Adding sequence "%s" to string set.', seq)
        string_set.append(seq)
    logging.info('Building index.')
    index = stempy.build_index(string_set)
    logging.info('Creating data object.')
    return stempy.Data(index)

def feq(x, y, eps=1e-4):
    return fabs(x - y < eps)



seq = "ACGTACACAC"
data = create_data(seq)
logging.info('Creating Markov model.')
mm, freqs = stempy.create_markov_model_order_3(data, 1.)
logging.info('Calculating likelihoods.')
lls = mm.calculate_likelihoods(data)
base_probs = map(exp, base_lls(lls[0]))
logging.info(', '.join('%.5f' % p for p in base_probs))
assert feq(stempy.W_mer_log_likelihood(lls[0], 0, 1), (seq.count('A')+1.)/(len(seq)+4.))
assert feq(base_probs[1], (seq.count('AC')+1.)/(seq.count('A')+4.))
assert feq(base_probs[2], (seq.count('ACG')+1.)/(seq.count('AC')-1.+4.)) # -1. because last 'AC' has no following character
assert feq(base_probs[3], (seq.count('ACGT')+1.)/(seq.count('ACG')+4.))
assert feq(base_probs[4], (seq.count('CGTA')+1.)/(seq.count('CGT')+4.))

# check the bg model from likelihoods
bg_model = stempy.create_bg_model_from_base_likelihoods(4, data, lls, freqs)

seq = "AAAC"
data = create_data(seq)
mm, freqs = stempy.create_markov_model_order_3(data, 1.)
assert feq(freqs.freq(0), .375)
assert feq(freqs.freq(1), .125)
assert feq(freqs.freq(2), .125)
assert feq(freqs.freq(3), .375)

seq = "AAACNNNNNNNTCTCTATACGCAGTACGG"
data = create_data(seq)
mm, freqs = stempy.create_markov_model_order_3(data, 1.)
lls = mm.calculate_likelihoods(data)
print ', '.join(map(str, lls))
for i, (x, y) in enumerate(pairs(lls[0])):
    assert x > y, '%d: %f <= %f' % (i, x, y)

