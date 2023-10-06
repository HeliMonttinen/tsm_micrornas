import os
import sys
import subprocess

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    Align a file with Pagan. At first create an initial
    alignment, then remove those nodes from the tree
    that are not present in the alignment.
    Then create ancestral sequence reconstruction.

    Alignments
    ==========

    alignment_dir: A path to a directory where the initial fasta
                   sequence is located.
    tree_dir: A path to a directory where the tree files are saved.
    subtree_align_dir: A path to a directory where subtrees
                       are saved.
    number: An Ensembl identifier.
    """

    alignment_dir = sys.argv[1]
    tree_dir = sys.argv[2]
    subtree_align_dir = sys.argv[3]
    number = sys.argv[4]

    identifier = number

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
             "{subtree_align_dir}{identifier}".format(
                 identifier=identifier,
                 subtree_align_dir=subtree_align_dir),
             "--ncbi",
             "--force-gap"])

    subprocess.call(
            ["/data/hmonttin/RNA/bash_scripts/pagan2",
             "-a",
             "{tree_dir}{identifier}.fas".format(
                 identifier=identifier,
                 tree_dir=tree_dir),
             "-r",
             "{tree_dir}{identifier}.tree".format(
                 identifier=identifier,
                 tree_dir=tree_dir),
             "-o",
             "{subtree_align_dir}{identifier}".format(
                 identifier=identifier,
                 subtree_align_dir=subtree_align_dir),
             "--guidetree"])

    subprocess.call(
            ["/data/hmonttin/RNA/bash_scripts/pagan2",
             "-a",
             "{tree_dir}{identifier}.fas".format(
                 identifier=identifier,
                 tree_dir=tree_dir),
             "-r",
             "{tree_dir}{identifier}.nhx_tree".format(
                 identifier=identifier,
                 tree_dir=tree_dir),
             "-o",
             "{subtree_align_dir}{identifier}_pagan".format(
                 identifier=identifier,
                 subtree_align_dir=subtree_align_dir),
             "--ancestors"])


if __name__ == "__main__":
    main()
