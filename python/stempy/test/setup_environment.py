#
# Copyright John Reid 2010, 2011
#

"""
Generic setup code shared between tests.
"""

import logging, os, sys

def init_test_env(f, level=logging.DEBUG):
    "Initialise test environment."
    logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")
    # logging.getLogger('').addHandler(logging.FileHandler('%s.log' % os.path.splitext(os.path.basename(file))[0]))
    logging.info('Command line: %s', ' '.join(sys.argv))

def append_to_path(d):
    logging.info('Appending to sys.path: %s', d)
    sys.path.append(d) # stempy
    
def prepend_to_path(d):
    logging.info('Prepending to sys.path: %s', d)
    sys.path.insert(0, d) # stempy
    
def update_path_for_stempy():
    d = os.path.dirname(__file__)
#     prepend_to_path(os.path.normpath(os.path.join(d, '..', '..', '..', '..', 'Infpy', 'python'))) # Infpy
#     prepend_to_path(os.path.normpath(os.path.join(d, '..', '..', '..', '..', 'PyICL', 'Python'))) # PyICL
#     prepend_to_path(os.path.normpath(os.path.join(d, '..', '..', '..', '..', 'Python', 'Cookbook', 'python'))) # cookbook
#     prepend_to_path(os.path.normpath(os.path.join(d, '..', '..'))) # stempy
    #print sys.path

def fasta_dir():
    return os.path.normpath(os.path.join(os.path.dirname(__file__), 'fasta'))

def is_debug_python():
    "@return: True iff is a debug python."
    return hasattr(sys, "gettotalrefcount") # only available in python debug build
