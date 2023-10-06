"""
Visualization of the secondary structures
with RNAplot. Colors the TSM source and target
sites on the structure.
"""

from collections import OrderedDict
import os
import sys
from pymongo import MongoClient
from subprocess import check_output
import shutil

from Bio import SeqIO

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    Arguments
    ==========

    target_dir: A path directory where plot.files are saved
    source_dir: Source directory of the .ps_files
    database: A database about TSM cases
    alignment_dir: A directory where sequence alignments
                   are located.

    Output
    =======
    _ss.ps files in the target directory
    """

    from RNA_alignments import indexes_in_original_seq

    target_dir = sys.argv[1]
    source_dir = sys.argv[2]
    database = sys.argv[3]
    alignment_dir = sys.argv[4]

    query = {"mismatches": {"$gt": 5},
             "TSM_length": {"$gt": 15},
             "quality": True}

    client = MongoClient()
    db = client[database]
    collection = db['TSMs']

    for document in collection.find(query):

        hit_dict = {}

        for key, value in document.items():

            hit_dict[key] = value

        identifier = hit_dict["identifier"]

        query_node = hit_dict["query"]
        ref_node = hit_dict["parent"]
        ancestor_node = hit_dict["ancestor"]

        bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                                    alignment_dir + identifier + '_pagan.fas',
                                    "fasta"))
        align_dict = OrderedDict()

        for seq in bio_align_dict:

            align_dict[seq] = str(bio_align_dict[seq].seq)
            if 'Homo' in seq:
                homo_sapiens_node = seq

        source_start = hit_dict["tsm_source_start"] + 1
        source_end = hit_dict["tsm_source_end"]

        target_start = hit_dict["tsm_target_start"] + 1
        target_end = hit_dict["tsm_target_end"]

        for node in [query_node, ref_node, homo_sapiens_node]:

            filename = source_dir + identifier + '_' + node + '.db'

            target_n = target_start

            if align_dict[node][target_start] == '-':
                target_n = target_start
                try:
                    while align_dict[node][target_n] == '-':
                        target_n += 1
                        print(target_n, align_dict[node][target_n])
                    target_n += 1
                except:
                    target_n = target_start

            new_indexes = indexes_in_original_seq(
                    [source_start, source_end, target_n, target_end],
                    align_dict[node])

            target_coloring = str(new_indexes[2]) + ' ' +\
                str(new_indexes[3]) + ' 12 1.00 0.60 0.67 omark'
            source_coloring = str(new_indexes[0]) + ' ' +\
                str(new_indexes[1]) + ' 12 0.72 0.69 0.47 omark'

            coloring = target_coloring + ' ' + source_coloring
            struct_file = filename

            check_output(
                    ["RNAplot",
                     "--layout-type=4",
                     "--pre",
                     "{coloring}".format(
                         coloring=coloring),
                     "{struct_file}".format(
                         struct_file=struct_file)])

            short_name_list = struct_file.rstrip('.db').split(
                    '/')[-1].split('_')
            short_name = '_'.join(short_name_list[1:]) + '_ss.ps'

            if node == query_node:
                tag = "qry"
            elif node == ref_node:
                tag = "ref"
            elif node == ancestor_node:
                tag = "anc"
            elif node == homo_sapiens_node:
                tag = "homo"

            target_file = target_dir + identifier + '_' +\
                str(hit_dict["position"]) + "_" + tag + "_" +\
                short_name.rstrip('_ss.ps') + '.ps'
            shutil.copyfile(short_name, target_file)

            if query_node == homo_sapiens_node:
                target_file = target_dir + identifier + '_' +\
                    str(hit_dict["position"]) + "_" + "homo" + "_" +\
                    short_name.rstrip('_ss.ps') + '.ps'
                shutil.copyfile(short_name, target_file)


if __name__ == "__main__":
    main()

