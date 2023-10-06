"""
This script renames tree nodes and make them to match
the sequence names. Finally it copies sequence file into
a new directory and writes tree_files with new node names.
"""

import os
import sys
import shutil

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    Rename the tree nodes.

    tree_file: A path to the tree file
    alignment_dir: A directory, where the original
    sequence files are located
    renamed_align_dir: A directory, where the sequence file is copies
    renamed_tree_dir: A directory, where the renamed tree is copied.
    """

    from common import return_file
    from tree_parse import rename_nodes_maf

    tree_file = sys.argv[1]
    alignment_dir = sys.argv[2]
    renamed_align_dir = sys.argv[3]
    renamed_tree_dir = sys.argv[4]

    for filename in return_file(alignment_dir):

        file_id = filename.rstrip('.fas').split('/')[-1]

        shutil.copyfile(filename.rstrip('.fas'),
                        renamed_align_dir + file_id)

        rename_nodes_maf(tree_file,
                         file_id,
                         renamed_tree_dir + file_id + '.tree',
                         alignment=alignment_dir,
                         align_outfile=renamed_align_dir + file_id)


if __name__ == "__main__":
        main()
