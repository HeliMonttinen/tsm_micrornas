"""
This tool is to identify distinct sequences in
a given set of alignments.
"""

from decimal import Decimal
import os
import sys
from ete3 import Tree

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    This script identifies sequences that are over
    2.5 times longer than human sequence, if the sequence contains
    more than 10 'N's and if the sequence identity between
    the human sequence and target sequence in non-gapped is
    50% or less, and the sequence identity between
    the human sequence and the target sequence in non-gapped
    regions of the target sequence is 60% or less.

    Arguments
    =========

    target_dir: A directory where cleaned sequence alignments
                are saved.

    alignment_dir: A directory where the original alignment files
                   are located.

    microRNA_inf: A file where information on all the human microRNA
                  loci is saved (data/210420_microRNAs_loc.txt).

    species_int: A species of interest to which all the other sequences
                 are compared.
    """

    from RNA_alignments import (parse_sequences_from_pairwise_mafft,
                                indexes_in_alignment)

    target_dir = sys.argv[1]
    alignment_dir = sys.argv[2]
    microRNA_inf = sys.argv[3]
    species_int = sys.argv[4]

    micro_info = {}

    with open(microRNA_inf, 'r') as f:

        for line in f:
            line_splitted = line.split()

            start = int(line_splitted[0].split(':')[1].split('-')[0])
            end = int(line_splitted[0].split('-')[1])

            ids = line_splitted[1]
            micro_info[ids] = [start, end]

    for subdir, dirs, files in os.walk(alignment_dir):
        for file in files:
            
            if 'pagan' in file or not file.endswith('.fas'):
                continue
            identifier = file.split('.fas')[0]

            tree = Tree(alignment_dir + identifier + '.nhx_tree')

            fasta_dict = parse_sequences_from_pairwise_mafft(
                    alignment_dir + identifier + '.fas')

            copy_dict = {}
            corr_start = None
            corr_end = None
            mistake = False

            for key in fasta_dict:
                spes_len = len(species_int.split('_'))
                if species_int in key and len(key.split('_')) == spes_len + 3:

                    start = int(key.split('_')[-2])
                    end = int(key.split('_')[-1])

                    if (start < micro_info[identifier][0]-50) or\
                            (end > micro_info[identifier][1]+50):
                        corr_start = (micro_info[identifier][0]) - 50 - start
                        corr_end = len(fasta_dict[key].replace('-', ''))\
                            - (end - (micro_info[identifier][-1] + 50))

                        if corr_end >= (len(fasta_dict[key].replace('-', ''))):
                            corr_list = [corr_start]
                        else:
                            corr_list = [corr_start, corr_end]

                        new_indexes = indexes_in_alignment(
                                corr_list,
                                fasta_dict[key])

                        if len(new_indexes) < 2:
                            mistake = True
                    elif (start > micro_info[identifier][0]-50) or\
                            (end < micro_info[identifier][-1]+50):
                        mistake = True

                    species_int_key = key
            if mistake is True:
                continue

            species_int_key = ""
            for key in fasta_dict:
                if species_int in key:
                    species_int_key = key
            if len(species_int_key) == 0:
                continue

            for key in fasta_dict:
                if len(fasta_dict[key].replace('-', '')) >\
                        len(fasta_dict[species_int_key].replace('-', ''))*1.3:
                    continue
                elif 'NNNNNNNNNN' in fasta_dict[key].replace('-', ''):
                    continue

                elif corr_start is not None:
                    
                    new_start = new_indexes[0]
                    new_end = 0
                    if len(corr_list) > 1:
                        new_end = new_indexes[1]
                    if len(corr_list) == 1:
                        corr_seq = fasta_dict[key][new_start:]
                    else:
                        corr_seq = fasta_dict[key][new_start:new_end]
                    if len(corr_seq.replace('-', '')) == 0:
                        continue
                    start_diff = len(
                            (fasta_dict[key][:new_start]).replace('-', ''))
                    if new_end != 0:
                        end_diff = len(
                                (fasta_dict[key][new_end:]).replace('-', ''))
                        new_end_ind = int(key.split('_')[-1]) - end_diff
                    else:
                        new_end_ind = int(key.split('_')[-1])
                    new_start_ind = int(key.split('_')[-2]) + start_diff
                    new_key = '_'.join(key.split('_')[:-2]) + '_' +\
                        str(new_start_ind) + '_' + str(new_end_ind)
                    copy_dict[new_key] = corr_seq
                    (tree&key).name = new_key
                    if species_int in new_key:
                        species_int_new = new_key

                else:
                    copy_dict[key] = fasta_dict[key]
                    if species_int in key:
                        species_int_new = key

            with open(target_dir + identifier, 'w') as f2:
                for item in copy_dict:

                    length = 0
                    count = 0
                    length2 = 0
                    count2 = 0

                    for a, b in zip(copy_dict[item],
                                    copy_dict[species_int_new]):

                        if a != '-' and b != '-':
                            length += 1
                            if b == a:
                                count += 1

                        if a != '-':
                            length2 += 1
                            if b == a:
                                count2 += 1

                    if length == 0 or Decimal(count/length) <= 0.5 or\
                            Decimal(count2/length2) < 0.6:
                        continue

                    f2.write('>' + item + '\n')
                    f2.write(copy_dict[item] + '\n')
            tree.write(outfile=target_dir + identifier + '.tree')


if __name__ == "__main__":
    main()
