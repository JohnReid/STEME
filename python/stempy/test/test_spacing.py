#
# Copyright John Reid 2012
#

"""
Test STEME motif spacing functionality.
"""

# do this first as otherwise have strange segmentation violation. Not sure of reason for this
import scipy.special

from setup_environment import init_test_env, logging
init_test_env(__file__, level=logging.INFO)

import pkg_resources
from optparse import OptionParser
from stempy.scan import load_occurrences_from_stream
from stempy.spacing import add_max_distance_option, count_all_pairs, spacing_idx

parser = OptionParser()
options = parser.get_default_values()
options.max_distance = 4


#
# Load the occurrences and associated sequence lengths,
# they will come sorted by position
#
logging.info('Loading occurrences')
occurrences, seq_infos, motifs = load_occurrences_from_stream(
    pkg_resources.resource_stream('stempy', 'test/spacing/steme-pwm-scan.out'),
    pkg_resources.resource_stream('stempy', 'test/spacing/steme-pwm-scan.seqs'),
)


#
# Iterate through the occurrences counting spacings
#
logging.info(
    'Examining spacings of up to %d b.p. between %d occurrences of %d motifs in %d sequences',
     options.max_distance, len(occurrences), len(motifs), len(seq_infos)
)
spacings = count_all_pairs(occurrences, seq_infos, ignore_close_to_end=True, options=options)

#
# Check we have the counts we expected
#
assert (spacings[('MOTIF-1', 'MOTIF-2')].sum() >= 0).all()
assert (spacings[('MOTIF-1', 'MOTIF-2')].sum() == 3).all()
assert spacings[('MOTIF-1', 'MOTIF-2')][spacing_idx(max_distance=4, distance=4, upstream=False, same_strand=True)] == 2
assert spacings[('MOTIF-1', 'MOTIF-2')][spacing_idx(max_distance=4, distance=4, upstream=True, same_strand=False)] == 1

assert (spacings[('MOTIF-2', 'MOTIF-1')].sum() >= 0).all()
assert (spacings[('MOTIF-2', 'MOTIF-1')].sum() == 3).all()
assert spacings[('MOTIF-2', 'MOTIF-1')][spacing_idx(max_distance=4, distance=4, upstream=True, same_strand=True)] == 2
assert spacings[('MOTIF-2', 'MOTIF-1')][spacing_idx(max_distance=4, distance=4, upstream=True, same_strand=False)] == 1

