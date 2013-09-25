#
# Copyright John Reid 2012
#

"""
Test find instances implementation.
"""

from setup_environment import init_test_env, logging
init_test_env(__file__, level=logging.INFO)

import os, sys, stempy, numpy


#
# Test parameters
#
Z_threshold = 1e-12
Z_threshold = .3
fasta = os.path.join(os.path.dirname(__file__), 'fasta', 'SOX2-t=20.fasta')
lambda_ = 0.006764069264069264
log_pwm = numpy.array(
    [
        [ -3.91087378,  -3.21835102,  -2.40776786,  -0.16258951],
        [-10.59673473,  -0.41557756,  -1.60941292,  -1.9660343 ],
        [ -1.04985069,  -2.40776786,  -2.99533239,  -0.67339553],
        [ -1.46966728,  -2.20714766,  -0.65397839,  -1.9660343 ],
        [ -1.30934073,  -2.30243512,  -1.27297639,  -1.04985069],
        [ -0.82102373,  -2.52551619,  -1.20398947,  -1.71475954],
        [ -3.50582491,  -2.12015522,  -1.71475954,  -0.40054025],
        [-10.59673473,  -0.79855214, -10.59673473,  -0.59789154],
        [ -2.12015522,  -1.9660343 ,  -2.52551619,  -0.41557756],
        [ -0.86754104, -10.59673473,  -0.89163714,  -1.77190979],
        [ -1.23788815,  -0.9942847 ,  -1.46966728,  -2.20714766],
        [ -1.42711219,  -1.46966728,  -1.20398947,  -1.46966728],
        [ -1.89705333,  -1.38629436,  -1.23788815,  -1.17120233],
        [ -2.40776786, -10.59673473,  -2.30243512,  -0.21079016]
    ]
)
W = len(log_pwm)





#
# Set up options
#
options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-find-instances-2')
options.bg_model_order = 4



#
# Load the sequences
#
input_sequences = stempy.SequenceSet(fasta.encode(sys.stdin.encoding or 'ascii'), options)



#
# Initialise the background
#
bg_manager = stempy.get_background_manager(input_sequences, input_sequences.mm, options)



#
# Create the model
#
model = input_sequences.create_model(bg_manager.get_bg_model(W), W)
model.bs.pssm.log_probs.values()[:] = log_pwm
model.bs.recalculate()
model.lambda_ = lambda_




#
# Create the instance finder and find the instances
#
instance_finder = stempy.FindInstances(input_sequences.data, model, Z_threshold)
instance_finder()
instance_finder.instances.sort()
logging.info('Found %d instances', len(instance_finder.instances))
