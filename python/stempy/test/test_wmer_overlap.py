#
# Copyright John Reid 2011
#

"""
Test STEME's performance on some small data sets of random sequences with planted sites.
"""

from setup_environment import init_test_env, logging
init_test_env(__file__, level=logging.INFO)

import stempy

evaluations1 = stempy.InstanceVec()
evaluations1.append(stempy.FindBestWMers.Evaluation(.5, 20, False))
evaluations1.append(stempy.FindBestWMers.Evaluation(.5,  0, False))
evaluations1.append(stempy.FindBestWMers.Evaluation(.5, 10, False))

evaluations2 = stempy.InstanceVec()
evaluations2.append(stempy.FindBestWMers.Evaluation(.5, 13, False))
evaluations2.append(stempy.FindBestWMers.Evaluation(.5, 23, False))
evaluations2.append(stempy.FindBestWMers.Evaluation(.5,  3, False))

stempy.sort_instances_by_position(evaluations1)
stempy.sort_instances_by_position(evaluations2)

for e in evaluations1:
    print e.global_pos    
print

for e in evaluations2:
    print e.global_pos
print
    
assert 3 + 6 + 6 == stempy.calculate_overlap(evaluations1, 6, evaluations2, 10)