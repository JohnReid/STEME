#
# Copyright John Reid 2010, 2011
#

"""
Test find best W-mers.
"""

from setup_environment import init_test_env, fasta_dir
init_test_env(__file__)

import stempy, logging, os, unittest
from cookbook.dicts import DictOf

stempy._dummy_fn()

def get_fasta_file(filename):
    return os.path.join(fasta_dir(), filename)



def check_w_mers_dont_overlap(W, w_mer_positions):
    w_mers_so_far = DictOf(list)
    for seq, pos in w_mer_positions:
        seen_so_far = w_mers_so_far[seq]
        for p in seen_so_far:
            if abs(p - pos) < W:
                raise ValueError('W-mers at positions %d and %d in sequence %d overlap' % (p, pos, seq))
        seen_so_far.append(pos)





class TestFindBestWMers(unittest.TestCase):
    
    def setUp(self):
        self.options = stempy.get_default_options()
        self.options.output_dir = os.path.join('output', 'test-find-best-w-mers')
        self.options.bg_model_order = 0
        
    def do_test(self, seed, fasta_file):
        "Do one test case."
                
        #
        # Load sequences and build index
        #
        algorithm = stempy.Algorithm(self.options)
        algorithm._initialise(fasta_file)
        data = algorithm.input_sequences.data

        #
        # look for best W-mers
        #
        W = len(seed)
        num_to_find = 3
        logging.info('Looking for %d best W-mers', num_to_find)
        best_w_mer_finder = stempy.create_best_w_mer_finder(
            data, 
            algorithm.create_model_of_input(W), 
            num_to_find
        )
        logging.info('Seeding model with %s', seed)
        if W != len(seed):
            raise ValueError('Seed must be same length as motif.')
        best_w_mer_finder.model.bs.seed(seed, True)
        best_w_mer_finder.model.set_lambda_for_sites(best_w_mer_finder.data.num_sequences)
        best_w_mer_finder()
        if not best_w_mer_finder.best_w_mers:
            raise RuntimeError('Did not find any W-mers')
        
        #
        # Log best W-mers
        #
        for _eval in best_w_mer_finder.best_w_mers:
            seq, offset = data.pos_localise(_eval.global_pos)
            strand = _eval.rev_comp and '-' or '+'
            logging.info(
                'Seed: %s; Site: %s; seq: % 2d; offset: % 4d; strand: %s; p(binding): %.2e; p(not binding): %.2e',
                seed, data.get_W_mer(W, _eval.global_pos), seq, offset, strand, _eval.Z, 1.-_eval.Z
            )

        #
        # check we at least found the seed...   
        #
        for _eval in best_w_mer_finder.best_w_mers:
            if data.get_W_mer(W, _eval.global_pos) == seed:
                break
        else:
            raise RuntimeError('Could not find seed in best W-mers')

        #
        # check we have no overlaps
        #
        localised_positions = [data.pos_localise(_eval.global_pos) for _eval in best_w_mer_finder.best_w_mers]
        check_w_mers_dont_overlap(W, localised_positions)
        assert not best_w_mer_finder.has_overlapping
    

    def test_find_best_W_mers(self):
        "Check we can find the seed in the best W-mers."    
        for seed, fasta_file in (
            ('CCAACCGG'        , get_fasta_file('find-best-w-mers-test.fa')),
            ('AGCGTGCGGCTTCTAC', get_fasta_file('T00759-tiny.fa')),
        ):
            self.do_test(seed, fasta_file)
            
          
    def test_specific_case(self):  
        """
        Check we get the correct W-mers.
        """
        fasta_file = os.path.normpath(get_fasta_file('T00759-tiny.fa'))
        seed = 'AAAACCCA'
        W = len(seed)
        num_sites = 4
        self.options.max_num_sites = num_sites
        self.options.min_num_sites = num_sites
        
        #
        # Load sequences and build index
        #
        algorithm = stempy.Algorithm(self.options)
        algorithm._initialise(fasta_file)
        data = algorithm.input_sequences.data

        model = algorithm.create_model_of_input(W)
        model.bs.seed(seed, True)
        model.set_lambda_for_sites(data.num_sequences)

        # look for best W-mers under model
        best_w_mer_finder = stempy.create_best_w_mer_finder(data, model, num_sites)
        best_w_mer_finder()
        if len(best_w_mer_finder.best_w_mers) < num_sites:
            if len(best_w_mer_finder.best_w_mers) != model.num_W_mers:
                raise ValueError('Did not find enough W-mers')
        
        # We want to get these W-mers
        # 
        # 2011-08-09 10:11:32,846 - INFO - Z=8.00e-02; pos=      313 +; AAAACCCA; AAAACCCA
        # 2011-08-09 10:11:32,846 - INFO - Z=4.37e-02; pos=      668 -; TGAGTTTT; AAAACTCA
        # 2011-08-09 10:11:32,846 - INFO - Z=1.37e-02; pos=      710 -; TGGTTTTC; GAAAACCA
        # 2011-08-09 10:11:32,846 - INFO - Z=1.37e-02; pos=      681 -; TGGTTCTT; AAGAACCA
        # 
        for wmer, (global_pos, rev_comp) in zip(best_w_mer_finder.best_w_mers, [(313, False), (668, True), (710, True), (681, True)]):
            if wmer.global_pos != global_pos and wmer.Z < model.calculate_Z(global_pos, rev_comp):
                raise ValueError('Got wrong W-mers')
    
    
    def test_we_get_all_W_mers_we_asked_for(self):
        """
        Check that we are not short-changed on the number of W-mers we asked for.
        """
        fasta_file = os.path.normpath(get_fasta_file('T00759-small.fa'))
        num_sites = [2, 4, 8, 16, 32]
        self.options.max_num_sites = max(num_sites)
        self.options.min_num_sites = min(num_sites)
        
        #
        # Load sequences and build index
        #
        algorithm = stempy.Algorithm(self.options)
        algorithm._initialise(fasta_file)
        data = algorithm.input_sequences.data

        for seed in (
            'GCTAGCTAGCGG',
            'ATGCAGAAAAATTAAG',
            'TTTAAAATACTTTAAA',
        ):
            # seed a model
            logging.info('Using seed %s', seed)
            W = len(seed)
            model = algorithm.create_model_of_input(W)
            model.bs.seed(seed, True)
            model.set_lambda_for_sites(data.num_sequences)
        
            for num_to_find in num_sites:
                # look for best W-mers under model
                best_w_mer_finder = stempy.create_best_w_mer_finder(data, model, num_to_find)
                best_w_mer_finder()
                if len(best_w_mer_finder.best_w_mers) < num_to_find:
                    if len(best_w_mer_finder.best_w_mers) != model.num_W_mers:
                        logging.warning('Found %d W-mers', len(best_w_mer_finder.best_w_mers))
                        logging.warning('%d W-mers available', model.num_W_mers)
                        logging.warning('Wanted %d W-mers', num_to_find)
                        raise ValueError('Did not find enough W-mers')


    def test_find_best_W_mers_2(self):
        """
        Run best W-mer finder and check we found the seeds we wanted.
        """
        self.options.min_num_sites = self.options.max_num_sites = num_to_find = 2
        
        # load data and create STEME object
        fasta_file = os.path.normpath(get_fasta_file('T00759-small.fa'))
        
        #
        # Load sequences and build index
        #
        algorithm = stempy.Algorithm(self.options)
        algorithm._initialise(fasta_file)
        data = algorithm.input_sequences.data

        for seed in (
            'ATGCAGAAAAATTAAG',
            'TTTAAAATACTTTAAA',
        ):
            # create and seed a model
            W = len(seed)
            model = algorithm.create_model_of_input(W)
            model.bs.seed(seed, True)
            model.set_lambda_for_sites(data.num_sequences)
        
            # look for best W-mers under model
            best_w_mer_finder = stempy.create_best_w_mer_finder(data, model, num_to_find)
            best_w_mer_finder()
            avg_Z = 0.
            for _eval in best_w_mer_finder.best_w_mers:
                logging.info(
                    'Seed: %s; Site: %s; p(binding): %.2e; p(not binding): %.2e',
                    seed, data.get_W_mer(W, _eval.global_pos), _eval.Z, 1.-_eval.Z
                )
                avg_Z += _eval.Z
            logging.info('Seed: %s; Average Z: %.6f', seed, avg_Z / len(best_w_mer_finder.best_w_mers))
            
            #
            # Check we found the seed
            #
            for _eval in best_w_mer_finder.best_w_mers:
                if data.get_W_mer(W, _eval.global_pos) == seed:
                    break
            else:
                raise RuntimeError('Could not find seed in best W-mers')
            
            #
            # Log the product of p-values
            #
            best_w_mer_finder.update_model(num_to_find, use_pseudo_counts=False)
            logging.info('Seed: %s; log PoP: %.6f', seed, algorithm.significance.log_product_p_values(model))
        




if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)
    
    # only run for release builds
    import sys
    _python_debug_build = hasattr(sys, "gettotalrefcount") # only available in python debug build
    if not _python_debug_build:
        #TestFindBestWMers('test_find_best_W_mers_2').debug()
        unittest.main()
    else:
        TestFindBestWMers('test_find_best_W_mers').debug()
        TestFindBestWMers('test_specific_case').debug()
        
    