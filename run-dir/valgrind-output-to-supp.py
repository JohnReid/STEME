#!/usr/bin/env python
#
# Copyright John Reid 2010
#

"""
Script to convert valgrind log to a suppressions file.
"""


import sys, uuid

counter = 0
_id = uuid.uuid1()

for l in sys.stdin:
    if l.startswith('=='):
        continue
    l = l.strip('\n')
    if '   <insert_a_suppression_name_here>' == l:
        l = '   supp-%s-%05d' % (_id, counter)
        counter += 1
    print l

