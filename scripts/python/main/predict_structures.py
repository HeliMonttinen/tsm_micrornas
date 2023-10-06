"""
Predict secondary with RNAfold structures structures
for the fasta files in a given directory.

"""

from collections import OrderedDict
import os
import sys
from Bio import SeqIO


dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    Go through files in a given directory and
    predict structures with RNAfold.

    Arguments
    =========
    alignment_dir : A directory where the fasta files are
                    located
    structure_dir: A directory where structure files will be saved.

    """

    from common import return_file
    from structure import (sequences_for_single_structural_analysis,
                           analyse_structure)

    alignment_dir = sys.argv[1]
    structure_dir = sys.argv[2]

    for filename in return_file(alignment_dir):
        if not filename.endswith('_pagan.fas'):
            continue
        bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                                       filename,
                                       "fasta"))

        align_dict = OrderedDict()

        struct_list = []
        align_id = filename.split('/')[-1].rstrip('_pagan.fas')
        for seq in bio_align_dict:
            align_dict[seq] = str(bio_align_dict[seq].seq)
            if not os.path.exists(
                    structure_dir + align_id + '_' + seq + '.fasta'):
                struct_list.append(seq)

        sequences_for_single_structural_analysis(structure_dir,
                                                 align_id,
                                                 align_dict,
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


if __name__ == "__main__":
    main()

