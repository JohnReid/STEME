#
# Copyright John Reid 2011
#

"""
Test find starts.
"""

"""
Code to test the STEME EM algorithm.
"""

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append('/home/john/Dev/MyProjects/Python/Cookbook')

import logging
logging.basicConfig(level=logging.INFO)

import stempy


#
# parse options and run algorithm
#
options, args = stempy.parse_options(stempy.add_options)
#options.epsilon = 0.
options.min_w = 8
options.max_w = 8
fasta = os.path.join(os.path.dirname(__file__), 'fasta', 'T00759-small.fa')
algorithm = stempy.Algorithm(options)
algorithm(fasta)


