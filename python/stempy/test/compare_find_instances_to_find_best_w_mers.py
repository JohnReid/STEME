#
# Copyright John Reid 2011, 2012, 2013
#

"""
Compare find instances to find best W-mers in terms of efficiency.
"""

from setup_environment import init_test_env, logging
init_test_env(__file__)


import stempy
from cookbook.timer import Timer

seed = 'ATAAAA'
fasta = '/home/john/Data/MO-MK-EB/unique_MK.fasta.masked'
#fasta = '/home/john/Data/MO-MK-EB/MO_MK_EB_shared.fasta.masked'

options = stempy.get_default_options()
W = options.min_w = options.max_w = len(seed)

# load the sequences
num_bases, seqs, ids, index = stempy.read_sequences(fasta, options)

# create the data object
with Timer(msg='build data'):
    data = stempy.Data(index, max_W=options.max_w)

# get the background
mm, freqs = stempy.create_markov_model_order_from_index_4(data.index, options.back_dist_prior)
freqs_with_pseudo_counts = freqs.add_pseudo_counts(options.back_dist_prior)
lls = mm.calculate_likelihoods(data)
bg_model = stempy.create_bg_model_from_base_likelihoods(W, data, lls, freqs_with_pseudo_counts)

# binding site model
bs_model = stempy.PssmBindingSiteModel(stempy.initialise_uniform_pssm(W, options.alphabet_size))
bs_model.seed(seed)

# whole model
model = stempy.Model(data, bs_model, bg_model, _lambda=0.)

Z_threshold = .3
with Timer(msg='find instances with Z>%f' % Z_threshold):
    instance_finder = stempy.FindInstances(data, model, Z_threshold)
    instance_finder()
    logging.info('Found %d instances', len(instance_finder.instances))


num_W_mers_to_find = 10000
with Timer(msg='find %d best W-mers' % num_W_mers_to_find):
    w_mer_finder = stempy.create_best_w_mer_finder(data, model, num_W_mers_to_find)
    w_mer_finder()
    logging.info('Found %d instances', len(w_mer_finder.best_w_mers))
    

def global_overlap(pos1, pos2, W):
    return abs(pos1 - pos2) < W

def get_non_overlapping(instances, W):
    instances.sort()
    instances.reverse()
    result = []
    for i in instances:
        for better in result:
            if global_overlap(i.global_pos, better.global_pos, W):
                break
        else:
            result.append(i)   
    return result

def get_non_overlapping_by_position(instances, W):
    instances.sort_by_position()
    result = []
    last = None
    for this in instances:
        if not last:
            last = this
        elif global_overlap(last.global_pos, this.global_pos, W):
            last = max(last, this)
        else:
            result.append(last)
            last = this
    if last:
        result.append(last)
    return result

def check_non_overlapping(instances, W):
    if instances.do_instances_overlap(W, already_sorted_by_position=False):
        raise RuntimeError('Instances overlap')

with Timer('get non-overlapping by position'): 
    non_overlapping = instance_finder.instances.get_best_non_overlapping(W, already_sorted_by_position=False)
    logging.info('Got %d non-overlapping W-mers by position', len(non_overlapping))
check_non_overlapping(non_overlapping, W)

def make_instance_vec(instances):
    result = stempy.InstanceVec()
    result.extend(instances)
    return result

non_overlapping2 = get_non_overlapping_by_position(instance_finder.instances, W)

