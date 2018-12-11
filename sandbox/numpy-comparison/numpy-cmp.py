#!/usr/bin/env python
"""
Script to test numpy's elementwise comparison.

It seems the test `None != a` fails on numpy 1.13.0 but succeeds on 1.12.1 (with a FutureWarning).
See: https://docs.scipy.org/doc/numpy-1.13.0/release.html
"""
import numpy as np
print('Numpy version: {}'.format(np.__version__))
a = np.zeros(4)
if a is not None:
  print('a is not None')
if None != a:
  print('None != a')
