import os
import sys
import subprocess

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    This script realigns the sequences and removes
    branches, which are not present in the alignment and\
            finally infers the ancestral sequences.

    Arguments
    ==========

    alignment_dir: A directory, where the original alignment
                    is located.
    tree_dir: A directory, where tree files are located

    target_align_dir: A directory, where resulting alignments
                      and trees are saved

    identifier: An Ensembl identifier

    ."""

    alignment_dir = sys.argv[1]
    tree_dir = sys.argv[2]
    target_align_dir = sys.argv[3]
    identifier = sys.argv[4]

    subprocess.call(
            ["/data/hmonttin/RNA/bash_scripts/pagan2",
             "-s",
             "{alignment_dir}{identifier}".format(
                 identifier=identifier,
                 alignment_dir=alignment_dir),
             "-t",
             "{tree_dir}{identifier}.tree".format(
                 identifier=identifier,
                 tree_dir=tree_dir),
             "-o",
             "{target_align_dir}{identifier}".format(
                 identifier=identifier,
                 target_align_dir=target_align_dir),
             "--guidetree"])

    subprocess.call(
            ["/data/hmonttin/RNA/bash_scripts/pagan2",
             "-a",
             "{alignment_dir}{identifier}.fas".format(
                 identifier=identifier,
                 alignment_dir=alignment_dir),
             "-r",
             "{tree_dir}{identifier}.nhx_tree".format(
                 identifier=identifier,
                 tree_dir=tree_dir),
             "-o",
             "{target_align_dir}{identifier}_pagan".format(
                 identifier=identifier,
                 target_align_dir=target_align_dir),
             "--ancestors"])


if __name__ == "__main__":
    main()
