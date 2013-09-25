#
# Copyright John Reid 2011
#

"""
Plots the speed at which MEME initialises its p-value tables for various values of max-sites.
"""

import pylab as P

speeds = {
    #    (2,2)    : 0,
    #    (2,4)    : 0,
    #    (2,8)    : 0.01,
    #    (2,16)   : 0,
    (2, 32): 0,
    (2, 64): 0.03,
    (2, 128): 0.23,
    (2, 256): 2.05,
    (2, 512): 17.32,
    (2, 1024): 137.16,
    (2, 2048): 1051.88,
}

max_sites = [maxsites for (minsites, maxsites) in speeds]
max_sites.sort()
durations = [speeds[(2, k)] for k in max_sites]
P.clf()
P.loglog(max_sites, durations, 's')
P.xlabel('max sites')
P.ylabel('seconds')
# P.show()
