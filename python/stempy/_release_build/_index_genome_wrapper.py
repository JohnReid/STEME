#
# Copyright John Reid 2012
#

import logging
logger = logging.getLogger(__name__)

from ._index_genome import *
logger.info(
    'Loaded genome version of STEME C++-python interface: max seqs=%d; max seq length=%d; max total length=%d',
    index_max_seqs, index_max_seq_length, index_max_text_length
)



def get_background_manager(sequences, mm, options):
    import stempy
    return stempy.MarkovBgModelManager(sequences, mm, options)
