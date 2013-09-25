#!/usr/bin/env python
#
# Copyright John Reid 2011
#

"""
Test find best W-mers speed.
"""

# import sys
# _python_debug_build = hasattr(sys, "gettotalrefcount") # only available in python debug build
# if _python_debug_build:
#     import logging, os
#     logging.basicConfig(level=logging.INFO)
# 
#     def append_to_path(dir):
#         logging.info('Appending to sys.path: %s', dir)
#         sys.path.append(dir) # stempy
# 
#     def update_path_for_stempy():
#         dir = os.path.dirname(__file__)
#         append_to_path(os.path.normpath(os.path.join(dir, '..', '..'))) # stempy
#         append_to_path(os.path.normpath(os.path.join(dir, '..', '..', '..', '..', 'Infpy', 'python'))) # Infpy
#         append_to_path(os.path.normpath(os.path.join(dir, '..', '..', '..', '..', 'PyICL', 'Python'))) # PyICL
# 
#     update_path_for_stempy()



#
# Set up the logging
#
import logging, os, sys, time
from cookbook.script_basics import setup_logging
setup_logging(__file__, level=logging.INFO)
import stempy


def get_fasta_file(filename):
    return os.path.join(os.path.dirname(__file__), '..', 'fasta', filename)

seed = 'TTTAAAATACTTTAAA'
num_to_find = 10000

options = stempy.get_default_options()
options.max_num_sites = options.min_num_sites = 10
options.min_w = options.max_w = W = len(seed)

#
# read in data
#
if 1 < len(sys.argv):
    fasta_file = sys.argv[1]
else:
    fasta_file = os.path.normpath(get_fasta_file('T00759-small.fa'))
algorithm = stempy.Algorithm(options)
algorithm.initialise(fasta_file.encode(sys.stdin.encoding or 'ascii'))
algorithm._initialise_p_value_tables()
data = algorithm.data


#
# Create model, set pseudo-counts and seed
#
logging.info('Creating model')
bg = algorithm._get_bg_model(W)
bs = stempy.PssmBindingSiteModel(stempy.initialise_uniform_pssm(W, algorithm.options.alphabet_size))
model = stempy.Model(data, bs, bg, _lambda=algorithm.options.lambda_)
model.bs.seed_pseudo_counts = options.starts_seed_pseudo_counts
logging.info('Seeding model with %s', seed)
if W != len(seed):
    raise ValueError('Seed must be same length as motif.')
model.bs.seed(seed, True)
model.set_lambda_for_sites(data.num_sequences)


for BestWMerFinder in [
    stempy.FindBestWMersMultiIndex, 
    # stempy.FindBestWMersSet, 
    stempy.FindBestWMersSortedVec
]:
    #
    # look for best W-mers under model
    #
    logging.info('Looking for best %d W-mers using ************ %s ************', num_to_find, BestWMerFinder)
    start = time.clock()
    best_w_mer_finder = BestWMerFinder(data, model, num_to_find)
    best_w_mer_finder()
    logging.info('Found best %d out of %d W-mers in %f seconds', len(best_w_mer_finder.best_w_mers), model.num_W_mers, time.clock() - start)
    
    
    #
    # Show some of the best W-mers
    #
    num_to_show = 5
    avg_Z = sum(_eval.Z for _eval in best_w_mer_finder.best_w_mers) / len(best_w_mer_finder.best_w_mers)
    logging.info('Seed: %s; Average Z: %.6f', seed, avg_Z)
    logging.info('Showing best %d W-mers', num_to_show)
    for _eval in best_w_mer_finder.best_w_mers[:num_to_show]:
        logging.info(
            'Seed: %s; Site: %s; p(binding): %.2e; p(not binding): %.2e',
            seed, data.get_W_mer(W, _eval.global_pos), _eval.Z, 1.-_eval.Z
        )
    logging.info('Showing worst %d of the best W-mers', num_to_show)
    for _eval in best_w_mer_finder.best_w_mers[-num_to_show:]:
        logging.info(
            'Seed: %s; Site: %s; p(binding): %.2e; p(not binding): %.2e',
            seed, data.get_W_mer(W, _eval.global_pos), _eval.Z, 1.-_eval.Z
        )
    
    


