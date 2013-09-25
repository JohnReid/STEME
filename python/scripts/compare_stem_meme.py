#
# Copyright John Reid 2009, 2010
#

"""
Code to compare STEME and MEME algorithms.
"""


# set up logging
import logging
from cookbook.script_basics import setup_logging
setup_logging(__file__, level=logging.INFO)

import os
import math
import stempy
import stempy.meme as meme
from stempy import parse_options, read_fasta, ensure_dir_exists, hamming_distance, reverse_complement


def dir_for_options(options, fasta):
    """
    Return a directory for the settings in options.
    """
    return os.path.join(
        'output',
        'comparison',
        'W=%02d-%02d-sites=%d-fasta=%s' % (
            options.min_w,
            options.max_w,
            options.max_num_sites,
            os.path.splitext(os.path.basename(fasta))[0]
        )
    )


def find_best_w_mers_for_seed(stem_algorithm, seed, options, log_level):
    "Find the best W-mers for the seed."
    num_to_find = options.max_num_sites
    logging.log(
        log_level, 'Finding best %d W-mers for seed %s', num_to_find, seed)
    _freqs, _freqs_with_pseudo_counts, model = stempy.Model(
        stem_algorithm.data,
        len(seed),
        options
    )
    logging.debug('Seed pseudo-count: %f', model.bs.seed_pseudo_counts)
    logging.debug('Model lambda: %f', model.lambda_)
    model.bs.seed(seed, True)
    best_w_mer_finder = stempy.create_best_w_mer_finder(
        stem_algorithm.data, model, num_to_find)
    best_w_mer_finder()
    best_w_mers = list(best_w_mer_finder.best_w_mers)[-num_to_find:]
    avg_Z = 0.
    avg_distance = 0.
    for eval_ in best_w_mers:
        w_mer = stem_algorithm.data.get_W_mer(len(seed), eval_.global_pos)
        if eval_.rev_comp:
            w_mer = reverse_complement(w_mer)
        distance = hamming_distance(seed, w_mer)
        logging.log(log_level, 'Best W-mer: %s; Z=%.5f, distance=%d',
                    w_mer, eval_.Z, distance)
        avg_Z += eval_.Z
        avg_distance += distance
    best_w_mer_finder.update_model(num_to_find, use_pseudo_counts=False)
    log_pop = model.log_product_of_pvalues
    avg_Z /= len(best_w_mers)
    logging.log(log_level, 'Seed: %s; log PoP: %.6f', seed, log_pop)
    logging.log(log_level, 'Seed: %s; Average Z=%.7f', seed, avg_Z)
    logging.log(log_level, 'Seed: %s; Average distance=%.2f',
                seed, avg_distance / len(best_w_mers))
    logging.debug('Frequencies: %s' %
                  ' '.join(map(str, (model.bg.freqs.freq(x) for x in xrange(4)))))
    return avg_Z, log_pop


def compare_meme_stem(options, fasta):
    """
    Runs MEME and STEME and compares results.
    """

    logging.info('Comparing STEME and MEME on %s', fasta)

    # run MEME
    meme_algorithm = meme.Algorithm(options)
    meme_algorithm(fasta)

    # run STEM
    stem_algorithm = stempy.Algorithm(options)
    stem_algorithm(fasta)

    # check if the starts were the same...
    if len(meme_algorithm.starts) != len(stem_algorithm.motifs[0].start_results):
        logging.warning(
            'STEME has different number of starts (%d) to MEME (%d)',
            len(stem_algorithm.motifs[0].start_results), len(meme_algorithm.starts))
    else:
        # check each start
        for meme_start_result, stem_start_result in zip(meme_algorithm.starts, stem_algorithm.motifs[0].start_results):
            meme_seed = meme_start_result.cons0
            meme_score = meme_start_result.sig
            stem_seed = stem_start_result.start.seed
            stem_score = math.exp(stem_start_result.log_E_value)

            if meme_score > stem_score:
                logging.warning(
                    'MEME has better score than STEME!!!!!!!!!!!!!!!!!!!! (%s) %f > (%s) %f', meme_seed, meme_score, stem_seed, stem_score)

            if meme_seed != stem_seed:

                # the start was different
                logging.info(
                    'STEME has different seed (%s) to MEME (%s)', stem_seed, meme_seed)

                # Examine the best W-mers for each start under STEM's model
                logging.info(
                    "Finding best W-mers for STEME's seed: %s", stem_seed)
                stem_avg_Z, stem_log_pop = find_best_w_mers_for_seed(
                    stem_algorithm, stem_seed, options, logging.INFO)
                logging.info(
                    "Finding best W-mers for MEME's seed: %s", meme_seed)
                meme_avg_Z, meme_log_pop = find_best_w_mers_for_seed(
                    stem_algorithm, meme_seed, options, logging.INFO)

                if meme_avg_Z > stem_avg_Z:
                    logging.warning(
                        "MEME's average Z better than STEME's!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

                if meme_log_pop < stem_log_pop:
                    logging.error(
                        "MEME's log product-of-p-values better than STEME's!!!!!!!!!!!!!!!!!!!!!!")

                # See if STEME converges to same motif when using MEME's start
#                model_score, model, EM, LLs = stem_algorithm.motifs[0]
#                model.bs.seed(meme_seed, True)
#                logo(N.exp(model.bs.pssm.log_probs.values()), 'test-em-0-seed', options.output_dir)
#                stempy.run_em_from_start(stem_algorithm.stem, model, meme_seed, options)
#                expected_sites = EM.do_iteration()
#                logo(N.exp(model.bs.pssm.log_probs.values()), 'test-em-1-after-EM', options.output_dir)

    return meme_algorithm, stem_algorithm


def num_sequences(fasta):
    "@return: The number of sequences in the fasta file."
    _num_bases, seqs, _ids = read_fasta(fasta)
    return len(seqs)


# parse options
def add_options(options):
    stempy.add_options(options)
    meme.add_options(options)
options, args = parse_options(add_options)
options.store_start_results = True

# for each fasta file
for fasta in [
    os.path.abspath(
        os.path.join(os.path.dirname(__file__), '../test/fasta/dm01r.fasta')),
    #os.path.abspath(os.path.join(os.path.dirname(__file__), '../test/fasta/simple-test-01.fa')),
    #os.path.abspath(os.path.join(os.path.dirname(__file__), '../test/fasta/simple-test-02.fa')),
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), '../test/fasta/T00759-tiny.fa')),
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), '../test/fasta/T00759-small.fa')),
]:
    num_seqs = num_sequences(fasta)

    # for each number of sites to use
    for num_sites in [
        2,
        5,
        10,
        20,
        50
    ]:
        # don't try if we don't have enough sequences
        if num_sites > num_seqs:
            continue

        options.min_num_sites = options.max_num_sites = num_sites

        # for each width
        for min_w, max_w in [
            #                ( 6,  6),
            #                ( 8,  8),
            #                (12, 12),
            (6, 16),
            #                (16, 16),
        ]:
            options.min_w = min_w
            options.max_w = max_w
            options.output_dir = dir_for_options(options, fasta)
            ensure_dir_exists(options.output_dir)
            logging.info(
                'Comparing MEME to STEM. Output dir=%s', options.output_dir)
            meme_algorithm, stem_algorithm = compare_meme_stem(options, fasta)
