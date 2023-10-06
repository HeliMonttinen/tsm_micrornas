"""
This script predicts secondary structures
for a sequence cut from a given alignment file.

"""


from Bio import SeqIO
from collections import OrderedDict
import os
import sys

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def microRNA_int(microRNA_info):
    """
    Creates a set of identifiers from
    a given .given csv file.

    Argument
    ========
    microRNA_info: A .csv file containing microRNA
                   identifiers in the first column.

    Returns
    =======

    microRNAs: A set of identifiers

    """

    microRNAs = set()
    with open(microRNA_info) as f:

        for line in f:

            line_splitted = line.rstrip().split('\t')
            microRNAs.add(line_splitted[0])

    return microRNAs


def cut_alignments(microRNAs_int,
                   alignment_dir,
                   output_dir,
                   microRNA_info):
    """
    Cuts alignments to match the human primicroRNA region.

    Arguments
    =========
    microRNAs_int: microRNA_set to bet studied
    alignment_dir: A directory where alignment files are saved
    output_dir: A directory where cut sequences are saved
    microRNA_info: A .csv file where microRNA coordinates are saved

    Returns
    =======

    align_dict: A dictionary about aligned sequences
    cut_start: A starting index for a cut sequence
    cut_end: A stopping index for a cut sequence
    microRNA: microRNA identifiers

    """

    from RNA_alignments import indexes_in_alignment
    from tree_parse import cut_align_seq_position

    species_start_stop = {}

    with open(microRNA_info, 'r') as f:

        for line in f:
            line_splitted = line.split()

            locus = line_splitted[-1]
            if (locus.rstrip()).endswith('_'):
                continue

            start = int((locus.split(':')[-1]).split('-')[-2])
            end = int((locus.split(':')[-1]).split('-')[-1])

            ids = line_splitted[0]
            species_start_stop[ids] = [start, end+1]

    for microRNA in microRNAs_int:

        pathname = alignment_dir + microRNA + '_pagan.fas'

        bio_align_dict = SeqIO.to_dict(SeqIO.parse(
            alignment_dir + microRNA + '_pagan.fas',
            "fasta"))

        align_dict = OrderedDict()

        for seq in bio_align_dict:

            align_dict[seq] = str(bio_align_dict[seq].seq)

        with open(pathname, 'r') as f:

            for line in f:

                if 'Homo' in line:
                    seq_int = line.rstrip().lstrip('>')

                    try:
                        head_start = int((seq_int).split('_')[-2])
                        head_end = int((seq_int).rstrip().split('_')[-1])
                    except:
                        head_end = len(align_dict[seq_int].replace('-', ''))
                        head_end = head_start + head_end

                    diff_start = species_start_stop[microRNA][0] - head_start-1
                    diff_end = len((align_dict[seq_int]).replace('-', '')) -\
                        (head_end - species_start_stop[microRNA][1])

        new_indexes = indexes_in_alignment(
                        [diff_start, diff_end], align_dict[seq_int])

        (minimum,
         maximum,
         cut_start,
         cut_end,
         cut_align_dict) = cut_align_seq_position(
                 align_dict,
                 new_indexes[0],
                 new_indexes[0],
                 new_indexes[1],
                 new_indexes[1])

        yield align_dict, cut_start, cut_end, microRNA


def predict_structs(cut_start,
                    cut_end,
                    align_dict,
                    structure_dir,
                    align_id):
    """
    Predict structures for cut squences with RNAfold.

    Arguments
    =========
    cut_start: A starting index of a cut alignment
    cut_end: A stopping index of a cut alignment
    align_dict: A dictionary about aligned sequences
    structure_dir: A directory where structures are saved
    align_id: An alignment identifier

    """

    from structure import (sequences_for_single_structural_analysis,
                           analyse_structure)

    struct_list = []
    cut_align_dict = {}
    for seq in align_dict:
        if '#' not in seq:
            cut_align_dict[seq] = align_dict[seq][cut_start:cut_end]

            struct_list.append(seq)

    sequences_for_single_structural_analysis(structure_dir,
                                             align_id,
                                             cut_align_dict,
                                             *struct_list)

    for struct in struct_list:

        try:
            analyse_structure(
                    structure_dir + align_id + '_' +
                    struct + ".fasta",
                    structure_dir + align_id + '_' +
                    struct + ".db",
                    "RNAfold")
        except:
            continue


def main():

    microRNAs_int = sys.argv[1]
    alignment_dir = sys.argv[2]
    output_dir = sys.argv[3]
    microRNA_info = sys.argv[4]

    microRNAs = microRNA_int(microRNAs_int)

    for align_dict, cut_start, cut_end, microRNA in cut_alignments(
            microRNAs,
            alignment_dir,
            output_dir,
            microRNA_info):

        predict_structs(
                cut_start,
                cut_end,
                align_dict,
                output_dir,
                microRNA)


if __name__ == "__main__":
    main()
