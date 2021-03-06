#
# Copyright John Reid 2010, 2011
#

"""
Test EM algorithm.
"""


#
# Trickery to find update path to import stempy from
#
from setup_environment import init_test_env, fasta_dir
init_test_env(__file__)

import stempy, os
from cookbook.named_tuple import namedtuple

Start = namedtuple('Start', 'seed num_sites score model best_w_mers')

options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-em-2')
seed = 'AAATTT'
W = len(seed)

# read the sequences and build STEME object from index
fasta = os.path.join(fasta_dir(), 'T00759-tiny.fa')
algorithm = stempy.Algorithm(options)
algorithm._initialise(fasta)
motif_finder = algorithm.create_motif_finder()

model = algorithm.create_model_of_input(W)
model.bs.seed(seed, True)
start = Start(seed=seed, num_sites=10, score=0., model=model, best_w_mers=stempy.InstanceVec())
em_result = motif_finder._run_em_from_start(start, return_EM=True)

# because we have symmetrical seed, should have same values for positive and negative strand Zs that exist
for n in xrange(algorithm.input_sequences.data.N):
    assert (
        not em_result.EM.get_Z(0).first
        and
        not em_result.EM.get_Z(0).second
    ) or (
        em_result.EM.get_Z(0).first.value == em_result.EM.get_Z(0).second.value
    )
