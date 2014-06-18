#
# Copyright John Reid 2014
#

"""
Test parse MEME minimal motif format
"""

from setup_environment import init_test_env, logging
init_test_env(__file__, level=logging.INFO)

import stempy.meme_parse as MP

filename = '../python/stempy/test/test-minimal-meme-format.txt'
text = open(filename).read()
parsed = MP.meme_format.parseString(text)
mmf = MP.extract_parsed(parsed)

mmf

