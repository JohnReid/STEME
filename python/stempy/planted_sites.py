#
# Copyright John Reid 2011
#

"""
Code to deal with generating sequences from a simple background model with planted sites in them.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Motif
from Bio.Alphabet import generic_dna
from itertools import imap
from collections import defaultdict
import numpy.random as R
import pyicl as P
import re
import logging


def idx_to_base(i):
    if 0 == i:
        return 'A'
    if 1 == i:
        return 'C'
    if 2 == i:
        return 'G'
    if 3 == i:
        return 'T'
    raise ValueError('Unknown index: %s' % i)


def seq_generator(avg_len):
    def gen_seq(s):
        length = R.poisson(avg_len)
        seq = Seq(
            "".join(imap(idx_to_base, R.randint(4, size=length))), generic_dna)
        record = SeqRecord(
            seq,
            id="Random-Seq-%02d" % s,
            description="A randomly generated sequence of length %d" % length)
        return record
    return gen_seq


def change_base(b, p_base_is_changed):
    if R.rand() < p_base_is_changed:
        return idx_to_base(R.randint(4))
    else:
        return b


def plant_sites(sequences, consensus, p_base_is_changed):
    logging.info(
        'Planting one binding site for consensus %s in each sequence', consensus)
    for s, seq in enumerate(sequences):
        site = Seq("".join(change_base(b, p_base_is_changed)
                   for b in consensus))
        pos = R.randint(len(seq) - len(consensus) + 1)
        mutable_seq = seq.seq.tomutable()
        rev_comp = bool(R.randint(2))
        site_info = 'site %s at position %3d (%s)' % (
            site, pos, rev_comp and "reverse complement" or "not reverse complement"
        )
        logging.info('Sequence %s: %s', seq.id, site_info)
        if rev_comp:
            site = site.reverse_complement()
        for i, b in enumerate(site):
            mutable_seq[i + pos] = b
        seq.seq = mutable_seq.toseq()
        seq.description += ': ' + site_info


_description_re = re.compile(
    'site ([ACGT]+) at position +(\d+) \((not )?reverse complement\)')


def parse_record_description(description):
    """
    Description will look like:
    
    "A randomly generated sequence of length 97: site AAAGGCTC at position 81 (reverse complement)"
    """
    match = _description_re.search(description)
    if None == match:
        raise ValueError('Could not parse description: %s', description)
    site = match.group(1)
    pos = int(match.group(2))
    rev_comp = None == match.group(3)
    return site, pos, rev_comp


def parse_fasta_for_sites(fasta_file):
    "Parse FASTA file to locate sites"
    true_sites = defaultdict(P.IntIntervalSet)
    logging.info('Parsing sites in %s', fasta_file)
    records = list(
        SeqIO.parse(open(fasta_file, "rU"), "fasta", alphabet=generic_dna))
    for record in records:
        site, pos, rev_comp = parse_record_description(record.description)
        logging.info(
            'Sequence: %s; site = %s; pos = %3d; rev_comp = %s', record.id, site, pos, rev_comp)
        true_sites[record.id].add(P.IntInterval(pos, pos + len(site)))
    return records, true_sites


def parse_meme_output_for_sites(meme_output):
    "Parse MEME-like output"
    logging.info('Parsing predictions from %s', meme_output)
    predicted_sites = defaultdict(P.IntIntervalSet)
    motifs = list(Motif.parse(open(meme_output), "MEME"))
    for motif in motifs:
        for instance in motif.instances:
            logging.info('Prediction: sequence = %s; site = %s; pos = %3d',
                         instance.sequence_name, instance, instance.start)
            predicted_sites[instance.sequence_name].add(
                P.IntInterval(instance.start, instance.start + len(instance)))
    return predicted_sites


def calculate_positives_and_negatives(records, true_sites, predicted_sites):
    "Gather true/false positives/negatives"
    TP = TN = FN = FP = 0
    for record in records:
        whole = P.IntInterval(0, len(record))
        true = true_sites[record.id]
        predicted = predicted_sites[record.id]
        TP += len(true.intersection(predicted))
        TN += len((true ^ whole).intersection(predicted ^ whole))
        FN += len(true.intersection(predicted ^ whole))
        FP += len((true ^ whole).intersection(predicted))
    return TP, TN, FN, FP
