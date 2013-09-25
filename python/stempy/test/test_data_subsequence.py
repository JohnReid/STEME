#
# Copyright John Reid 2010, 2011
#

"""
Test subsequence in data.
"""

from setup_environment import init_test_env, fasta_dir
init_test_env(__file__)

import stempy, logging, os, unittest

stempy._dummy_fn()

def get_fasta_file(filename):
    return os.path.join(fasta_dir(), filename)





class TestDataSubsequence(unittest.TestCase):
    
    def setUp(self):
        self.options = stempy.get_default_options()

    def test_data_subsequence(self):
        # read in data
        fasta_file = os.path.normpath(get_fasta_file('T00759-small.fa'))
        _num_bases, _seqs, _ids, index = stempy.read_sequences(fasta_file, self.options)
        data = stempy.Data(index)
        assert 'AGAGCG' == data.subsequence(2, 3, 6), 'AGAGCG != %s' % data.subsequence(2, 3, 6)


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)
    #TestDataSubsequence('test_data_subsequence').debug()
    unittest.main()
    