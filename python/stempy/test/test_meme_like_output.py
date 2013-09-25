#
# Copyright John Reid 2011
#

"""
Test STEME can produce MEME-like output.
"""

from setup_environment import init_test_env, logging
init_test_env(__file__, level=logging.INFO)

import os, stempy


#
# First run STEME
#
options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-meme-like-output')
options.min_w = options.max_w = 8
options.meme_like_output = 'meme.out'
algorithm = stempy.Algorithm(options)
fasta = os.path.join(os.path.dirname(__file__), 'fasta', 'T00759-tiny.fa')
algorithm(fasta)
logging.info('Showing MEME output from %s', algorithm.meme_like_output_file)
os.system('cat %s' % algorithm.meme_like_output_file)


#
# Test BioPython parser
#
from Bio import Motif
motifs = list(Motif.parse(open(algorithm.meme_like_output_file), "MEME"))


#
# Doesn't quite work with pycogent yet. Pycogent expects a summary section
# that contains sites in all the sequences 
#
#from cogent import LoadSeqs
#from cogent.parse.meme import MemeParser
#results = MemeParser(open(algorithm.meme_like_output_file, 'U'))
#seqs = LoadSeqs(fasta, aligned=False)
#results.Alignment = seqs
#for motif in results.Motifs:
#    module = motif.Modules[0]
#    print module.ID, module.Evalue, len(module.NamedSeqs)
#
#module_1 = results.Motifs[0].Modules[0]
#module_1.ID
#module_1.ConsensusSequence
#module_1.majorityConsensus(transform=str)
#module_1.IUPACConsensus()
#iupac = module_1.IUPACConsensus()
#majority = module_1.majorityConsensus()
#uncertainty = module_1.uncertainties()
#for i,m,u in zip(iupac,majority,uncertainty):
#    print i,m,u