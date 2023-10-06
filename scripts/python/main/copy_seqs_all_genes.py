"""
Extracts the full sequence from genomes based on the
given genomic coordinates
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
    the human sequence and the target sequence in non-gapped
    regions of the target sequence is 30% or less.

    Arguments
    =========

    target_dir: A directory where cleaned sequence alignments
                are saved.

    alignment_dir: A directory where the original alignment files
                   are located.

    identifier: A gene chunk identifier.
    """

    from RNA_alignments import (format_sequences_for_fpa,
                                parse_sequences_from_pairwise_mafft)

    target_dir = sys.argv[1]
    alignment_dir = sys.argv[2]
    identifier = sys.argv[3]

    tree = Tree(alignment_dir + identifier + '.tree')

    copy_dict = {}

    fasta_dict = parse_sequences_from_pairwise_mafft(
            alignment_dir + identifier + '.fas')

    found = False
    for item in fasta_dict:
        if "Homo" in item:
            Homo_sapiens = item
            copy_dict[Homo_sapiens] = fasta_dict[item]
            found = True
            break
    if found is False:
        return

    for key in fasta_dict:

        seq1, seq2 = format_sequences_for_fpa(
                fasta_dict[key],
                fasta_dict[Homo_sapiens])

        count = 0
        for a, b in zip(seq1, seq2):
            count += 1
            if a != '-' and b != '-':
                break
        if count > len(
                fasta_dict[Homo_sapiens].replace(
                    '-', '').replace('N', ''))*0.3:
            print(identifier, key, 'different lengths')
            continue

        count = 0
        for a, b in zip(seq1[:-1], seq2[:-1]):
            count += 1
            if a != '-' and b != '-':
                break
        if count > len(
                fasta_dict[Homo_sapiens].replace(
                    '-', '').replace('N', ''))*0.3:
            print(identifier, key, 'different lengths')
            continue

        elif 50*'N' in fasta_dict[key].replace('-', ''):
            print(identifier, "too many Ns")
            continue
        copy_dict[key] = fasta_dict[key]

    with open(target_dir + identifier, 'w') as f2:
        for item in copy_dict:

            length = 0
            count2 = 0
            for a, b in zip(copy_dict[item],
                            copy_dict[Homo_sapiens]):
                if a != '-':
                    length += 1
                    if a == b and a != 'N':
                        count2 += 1

            if length == 0:
                continue

            if Decimal(count2/length) < 0.30:

                continue

            f2.write('>' + item + '\n')
            f2.write(copy_dict[item] + '\n')
    tree.write(outfile=target_dir + identifier + '.tree')


if __name__ == "__main__":
    main()
