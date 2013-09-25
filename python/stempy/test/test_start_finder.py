#
# Copyright John Reid 2010, 2011
#

"""
Test start finding and associated code.
"""

from setup_environment import is_debug_python, init_test_env
init_test_env(__file__)

import stempy, logging, unittest

class TestStartFinder(unittest.TestCase):

    def test_edit_distances(self):  
        """
        Check our edit distances work.
        """
        assert  2 == stempy.edit_distance('AACC', 'GGCC')
        assert  2 == stempy.edit_distance('AACC', 'GAAC')
        assert  1 == stempy.shifted_edit_distance('AACC', 'GAAC')
        assert 16 == stempy.edit_distance('ACGTACGTACGTACGT', 'CGTACGTACGTACGTA')
        assert  1 == stempy.shifted_edit_distance('ACGTACGTACGTACGT', 'CGTACGTACGTACGTA')
        assert  4 == stempy.edit_distance('AACC', 'GGTT')
        assert  0 == stempy.rev_comp_shifted_edit_distance('AACC', 'GGTT')
        assert  0 == stempy.rev_comp_shifted_edit_distance('TTTC', 'GAAA')
        assert  1 == stempy.rev_comp_shifted_edit_distance('ATTTC', 'GAAAG')



if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)
    
    if is_debug_python():
        TestStartFinder('test_edit_distances').debug()
    else:
        unittest.main()
        
