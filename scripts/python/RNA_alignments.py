"""
Tools for aligning RNA sequences, identifying loops in the
structure and getting sequences.

Author: Heli Mönttinen (OrcID: 0000-0003-2461-0690)
"""

from collections import defaultdict
from decimal import Decimal
import subprocess
import numpy
import fpa_ext62

from common import split_multiple_fasta_file


def unusual_iupac_char(seq1, seq2):
    """
    Change chars in sequence pairs to standard one,
    if iupac char ranges overlap.

    Arguments
    ---------
    :seq1: sequence1
    :seq2: sequence2

    Returns
    -------
    :seq1: seq1 with standard char
    :seq2: se2 with standard char
    """

    iupac = defaultdict(list)
    iupac['R'].extend(['A', 'G'])
    iupac['Y'].extend(['C', 'T'])
    iupac['S'].extend(['G', 'C'])
    iupac['W'].extend(['A', 'T'])
    iupac['K'].extend(['G', 'T'])
    iupac['M'].extend(['A', 'C'])
    iupac['B'].extend(['C', 'G', 'T'])
    iupac['D'].extend(['A', 'G', 'T'])
    iupac['H'].extend(['A', 'C', 'T'])
    iupac['V'].extend(['A', 'C', 'G'])
    iupac['N'].extend(['A', 'C', 'G', 'T'])

    nseq1 = seq1.upper()
    nseq2 = seq2.upper()

    seq1 = ""
    seq2 = ""

    for a, b in zip(nseq1, nseq2):

        if a not in iupac and b not in iupac:

            seq1 = seq1 + a
            seq2 = seq2 + b
        elif a in iupac:
            a_set = set(iupac[a])
            if b in iupac:
                b_set = set(iupac[b])
            else:
                b_set = set([b])

            intersection = a_set.intersection(b_set)
            if len(intersection) > 0:
                common_char = list(intersection)[0]
                seq1 = seq1 + common_char
                seq2 = seq2 + common_char
            else:
                seq1 = seq1 + a
                seq2 = seq2 + b

        elif b in iupac:
            b_set = set(iupac[b])
            a_set = set([a])

            intersection = a_set.intersection(b_set)
            if len(intersection) > 0:
                common_char = list(intersection)[0]
                seq1 = seq1 + common_char
                seq2 = seq2 + common_char
            else:
                seq1 = seq1 + a
                seq2 = seq2 + b

    return seq1, seq2


def format_sequences_for_fpa(seq1, seq2):
    """
    Formats two sequences from multiple sequence alignment
    to be compatible with fpa. Changes "U" to "T"s and if both
    sequences have '-' in the same position, removes them.

    Requires:
    seq1 and seq2: Two sequences to compare from the multiple
                   sequence alignment.

    Returns:

    seq1, seq2 : Now in the right format

    """

    ns1 = seq1.lower().replace("u", "t")
    ns2 = seq2.lower().replace("u", "t")
    chr_list = set(['-', 'a', 't', 'c', 'g', 'r', 'y', 's',
                    'w', 'k', 'm', 'b', 'd', 'h', 'v', 'n'])

    seq1 = ""
    seq2 = ""

    for a, b in zip(ns1, ns2):
        if a == '-' and b == '-':
            continue
        else:
            if a not in chr_list:
                a = 'n'
            if b not in chr_list:
                b = 'n'
            a = a.upper()
            b = b.upper()
            seq1 = seq1 + a
            seq2 = seq2 + b

    return seq1, seq2


def parse_sequences_from_pairwise_mafft(stdout):
    """
    Takes a fasta-formatted output and parses
    the sequences into a dictionary format.
    """

    fasta_dict = {}

    for fasta in split_multiple_fasta_file(stdout):
        fasta_parts = fasta.split('\n')
        fasta_dict[fasta_parts[0].lstrip('>')] = fasta_parts[1]

    return fasta_dict


def run_fpa(seq1, seq2,
            perfect_copy=True,
            min_length=15,
            scan_window_limit=5):
    """
    runs fpa.
    Prints significat template_switch cases.

    Requires:

    seq1: First sequence to compare
    seq2: Second seqeunce to compare

    """
    FPA = fpa_ext62.FPA2()

    FPA.set_int("min_length", min_length)
    FPA.set_bool("clean_rna", True)
    FPA.set_int("scan_window_limit", scan_window_limit)
    FPA.set_int("ref_flank", 20)
    FPA.set_int("switch_flank", 20)
    FPA.set_int("scan_flank", 100)

    FPA.set_bool("maximise_23_score", True)
    FPA.set_bool('iupac', True)

    try:
        hits = FPA.scan_two(seq1, seq2, False)
    except:
        return '_', []

    for h in hits:

        hit_result = h.split(",")
        for i in numpy.arange(start=1, stop=len(hit_result)):
            hit_result[i] = Decimal(hit_result[i])

        yield h, hit_result


def fpa_print(raw_result, hit_result, id1, id2, seq1, seq2,
              loops=None, loop_seq=None, mode=None,
              min_length=15, scan_window_limit=5):
    """
    Formats the output of the fpa. Selects only interesting cases.
    Adds information on the RNA.
    Returns seqparate String streams for all_alignments as well
    as for those ones that have loops (mode='loops').

    Arguments
    =========
    hits1: The output from the run_fpa (seq1 as query).
    id1: identifier of the sequence.
    seq1: The first sequence
    loops1: The indexes of loops
    loop_seq1: Loop sequences.

    Those paramteres with '2' are the same as above, except they are
    characters of the second sequence.

    Returns
    =======

    fpa output, fpa output in loop mode (if mode is 'loop')

    """

    all_alignments = io.StringIO()
    proc = None
    """
    if hit_result[11] >= 0.85 and hit_result[10] >= 0.80 and\
           hit_result[12] >= 0.80 and  hit_result[11] >= 0.85 and\
           hit_result[13] > hit_result[15]:
    """
    if hit_result[11] >= 0.8 and\
            hit_result[13] > hit_result[15]:

        title = "#Sequences  " + id1 + ',' + id2 + '\n'
        all_alignments.write(title)

        all_alignments.write(raw_result + '\n')

        try:
            proc = subprocess.check_output(
                    ["python",
                     "-c",
                     "import fpa_ext62; fpa=fpa_ext62.FPA2();\
                     fpa.set_int('min_length', {min_length});\
                     fpa.set_bool('clean_rna', True); \
                     fpa.set_int('scan_window_limit', {scan_window_limit});\
                     fpa.set_int('ref_flank', 20);\
                     fpa.set_int('switch_flank', 20);\
                     fpa.set_int('scan_flank', 100);\
                     fpa.set_bool('maximise_23_score', True);\
                     fpa.set_bool('iupac', True);\
                     fpa.print_hit('{seq1}','{seq2}',\
                     '{raw_result}', True)".format(
                         seq2=seq2, seq1=seq1, min_length=min_length,
                         scan_window_limit=scan_window_limit,
                         raw_result=raw_result)])
        except:
            return None

        all_alignments.write(str(proc, 'utf-8'))

    if proc is not None:
        return all_alignments

    return None


def indexes_in_alignment(seq_indexes, aligned_seq):
    """
    Takes a list of indexes in the original sequence.
    These indexes are searched from the alignment.
    The function returns a list of corresponding indexes
    in the alignment.

    Parameters
    ----------

    :seq_indexes: A list of indexes from the original sequence
    :aligned_seq: The aligned sequence without identifier


    Result
    -------

    :alignment_indexes: A list of indexes in the alignment

    """

    alignment_indexes = []

    seq_len = len(str(aligned_seq).replace('-', '').replace('+', ''))

    for index in seq_indexes:

        int_index = index
        if index == seq_len:
            int_index = index - 1

        position_count = 0
        new_index = 0
        position_count = -1

        for position in aligned_seq:

            if position == '-' or position == '+':
                new_index += 1

            else:
                position_count += 1
                if position_count == int_index:
                    alignment_indexes.append(new_index)
                    break
                new_index += 1

    return alignment_indexes


def indexes_in_original_seq(aligned_indexes, aligned_seq):
    """
    Identifies the position of the loop in the original sequence.
    Alignment should be in one line. without header.

    Requires
    ========

    aligned_seq: Aligned sequence
    aligned_loop_indexes: Loop indexes in the alignment

    Returns
    =======

    Loop_indexes_original: Loop indexes in the original sequence

    """

    original_indexes = []

    for index in aligned_indexes:

        position_count = 0
        new_index = index

        for position in aligned_seq:
            if position == '-':
                new_index -= 1
                position_count += 1

            else:
                position_count += 1

            if position_count == index:
                original_indexes.append(new_index)
                break

    return original_indexes
