#
# Copyright John Reid 2010
#

"""
Code to parse MEME output like:

start_point: score 259.253166 cons GCAGAGGC
component  1: lambda= 0.093431 ps=     1370 cons0=GCAGAGGC score=259.25
Simplified        A  22525222
pos.-specific     C  25222225
probability       G  52252552
matrix            T  22222222


(start)   8  128.0 GCAGAGGC --> AGGGAGGA w   8 nsites    2 sig          156877.1709
"""


import sys, re
from cookbook.named_tuple import namedtuple

Start = namedtuple('Start', 'lambda_ ps score w0 nsites0 cons0 cons w nsites sig')

component_re = re.compile('component  [0-9]+: lambda=([ 0-9.]+)ps=([ 0-9.]+)cons0=[ACGT]+ score=([ 0-9.]+)')
start_re = re.compile('\(start\)   ([0-9]+) ([ 0-9.]+)([ ACGT]+)-->([ ACGT]+)w([ 0-9]+)nsites([ 0-9.]+)sig([ 0-9.]+)')

# for each line in input
for line in sys.stdin:
    
    if line.startswith('component'):
        match = component_re.match(line.strip())
        if not match:
            raise ValueError('Could not match "%s" to regular expression' % line)
        #print match.groups()
        lambda_ = float(match.group(1))
        ps = float(match.group(2))
        score = float(match.group(3))
        #print lambda_, ps, score
    
    elif line.startswith('(start)'):
        match = start_re.match(line.strip())
        if not match:
            raise ValueError('Could not match "%s" to regular expression' % line)
        #print match.groups()
        w0 = int(match.group(1))
        nsites0 = float(match.group(2))
        cons0 = match.group(3).strip()
        cons = match.group(4).strip()
        w = int(match.group(5))
        nsites = float(match.group(6))
        sig = float(match.group(7))
        #print w0, nsites0, cons0, cons, w, nsites, sig
        
        start = Start(lambda_=lambda_, ps=ps, score=score, w0=w0, nsites0=nsites0, cons0=cons0, cons=cons, w=w, nsites=nsites, sig=sig)
        print start
        

