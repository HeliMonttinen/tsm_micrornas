"""
This tool is for running fpa for all
parent-child node pairs within a given phylogenetic tree.
The identifier links alignment and tree files,
e.g. identifier.fas and identifier.anctree.

Before running copy fpa_ext6 file into a directory, in which
you are running your script.

"""

import os
import sys

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    This script is designed for going through a rooted tree
    and to run a FPA tool between each parent-child pair.
    As a result it writes a fpa result file (xxx_fpa) for a given tree.
    Each fpa result has a heading indicating the nodes compared:

    #Sequences  parental_node,child_node

    The alignment file has to be in a format identifier_pagan.fas and
    it has to contain alignments for internal nodes. The sequence identifiers
    have to correspond to node names in the tree.

    And the tree in a format identifier_pagan.anctree

    Arguments
    =========

    identifier: An identifier for the studied cluster.
    fpa_output_dir: A directory where the fpa results are written
    alignment_dir: A directory where the alignment files
                   (identifier.fasta) are located
    tree_dir: A directory where tree files are located
              (identifier_pagan.anctree).


    Output
    ======

    A file in fpa_output dir.
    """

    from RNA_alignments import (format_sequences_for_fpa,
                                run_fpa,
                                parse_sequences_from_pairwise_mafft,
                                fpa_print)

    from tree_parse import (read_tree,
                            go_through_branches)

    identifier = sys.argv[1]
    fpa_output_dir = sys.argv[2]
    alignment_dir = sys.argv[3]
    tree_dir = sys.argv[4]

    identifier = identifier.rstrip('_pagan.fas').rstrip('.fas')
  
    if not os.path.exists(alignment_dir + identifier + '_pagan.fas'):
        return

    fasta_dict = parse_sequences_from_pairwise_mafft(
        alignment_dir + identifier + '_pagan.fas')
    tree = read_tree(
        tree_dir + identifier + '_pagan.anctree')

    for anc_seqid, sample_seq1, sample_seq2 in go_through_branches(tree):
        anc_seq1, seq1 = format_sequences_for_fpa(
            fasta_dict[anc_seqid],
            fasta_dict[sample_seq1])

        anc_seq2, seq2 = format_sequences_for_fpa(
            fasta_dict[anc_seqid],
            fasta_dict[sample_seq2])
        for h, hit in run_fpa(anc_seq1, seq1):

            output1 = fpa_print(h, hit, anc_seqid,
                                sample_seq1, anc_seq1, seq1)

            if output1 is not None:
                o = output1.getvalue()
                with open(fpa_output_dir + identifier, 'a+') as f:
                    f.write(o)
            else:
                continue

        for h, hit in run_fpa(anc_seq2, seq2):
            output2 = fpa_print(h, hit, anc_seqid,
                                sample_seq2, anc_seq2, seq2)

            if output2 is not None:
                o = output2.getvalue()
                with open(fpa_output_dir + identifier, 'a+') as f:
                    f.write(o)
            else:
                continue


if __name__ == "__main__":
    main()
