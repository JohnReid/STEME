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
Z_threshold = .73 # this threshold works for both genome and standard build
#Z_threshold = .74985 # this threshold only works for standard build
W = 11
fasta = os.path.join(os.path.dirname(__file__), 'fasta', 'random-seqs-10-100.fasta')
log_pwm = numpy.array([
    [-0.2186892 , -4.02535169, -1.82812711, -4.02535169],
    [-0.05505978, -4.02535169, -4.02535169, -4.02535169],
    [-1.46040233, -0.65805586, -1.46040233, -4.02535169],
    [-4.02535169, -0.05505978, -4.02535169, -4.02535169],
    [-4.02535169, -4.02535169, -4.02535169, -0.05505978],
    [-1.82812711, -2.41591378, -0.98082925, -0.98082925],
    [-0.52884413, -4.02535169, -1.19213835, -2.41591378],
    [-1.82812711, -4.02535169, -0.41443378, -1.82812711],
    [-1.82812711, -1.82812711, -0.41443378, -4.02535169],
    [-0.41443378, -1.19213835, -4.02535169, -4.02535169],
    [-1.82812711, -2.41591378, -0.65805586, -1.46040233],
])




#
# Set up options
#
options = stempy.get_default_options()
options.output_dir = os.path.join('output', 'test-find-instances')



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




#
# Create the instance finder and find the instances
#
instance_finder = stempy.FindInstances(input_sequences.data, model, Z_threshold)
instance_finder()
instance_finder.instances.sort()
logging.info('Found %d instances', len(instance_finder.instances))
assert 13 <= len(instance_finder.instances)
# at least 13 instances in sequences
#2012-06-16 11:32:58,686 - INFO - seq=    5; pos=    67; strand=+; W-mer=AACCTCGAGAG; Z=0.749857
#2012-06-16 11:32:58,686 - INFO - seq=    0; pos=    48; strand=+; W-mer=AACCTAAGAAA; Z=0.814953
#2012-06-16 11:32:58,686 - INFO - seq=    3; pos=    51; strand=+; W-mer=AAACTGTGGCT; Z=0.819370
#2012-06-16 11:32:58,686 - INFO - seq=    5; pos=    79; strand=+; W-mer=AAGCTAAAGAG; Z=0.827948
#2012-06-16 11:32:58,687 - INFO - seq=    3; pos=    36; strand=-; W-mer=AAGCTTATCAG; Z=0.862206
#2012-06-16 11:32:58,687 - INFO - seq=    5; pos=    97; strand=-; W-mer=GAACTGGGGAT; Z=0.912242
#2012-06-16 11:32:58,687 - INFO - seq=    2; pos=    47; strand=+; W-mer=AAACTTGGGAA; Z=0.919355
#2012-06-16 11:32:58,687 - INFO - seq=    1; pos=     6; strand=+; W-mer=AACCTTAGACG; Z=0.963969
#2012-06-16 11:32:58,687 - INFO - seq=    6; pos=    46; strand=-; W-mer=AAGCTGGGGAC; Z=0.968255
#2012-06-16 11:32:58,687 - INFO - seq=    9; pos=    73; strand=-; W-mer=GACCTGATGAG; Z=0.968813
#2012-06-16 11:32:58,687 - INFO - seq=    5; pos=    16; strand=-; W-mer=AACCTGAGCCG; Z=0.974733
#2012-06-16 11:32:58,687 - INFO - seq=    6; pos=    73; strand=+; W-mer=AACCTTAGGCG; Z=0.984093
#2012-06-16 11:32:58,687 - INFO - seq=    3; pos=    10; strand=+; W-mer=AACCTTAGGAT; Z=0.984245



#
# Print the instances
#
for instance in instance_finder.instances:
    seq, pos = input_sequences.data.pos_localise(instance.global_pos)
    W_mer = input_sequences.data.get_W_mer(W, instance.global_pos)
    if instance.rev_comp:
        W_mer = stempy.reverse_complement(W_mer)
    logging.info('seq=%5d; pos=%6d; strand=%s; W-mer=%s; Z=%4f', seq, pos, instance.rev_comp and '-' or '+', W_mer, instance.Z)

