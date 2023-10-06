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
    species_int: The species, which microRNAs are studied
    splitted_len: The length of a splitted label of a species,
                  which microRNAs are studied. The character
                  used for splitting is '_'.

    Output
    =======
    _ss.ps files in the target directory
    """

    from RNA_alignments import indexes_in_original_seq

    target_dir = sys.argv[1]
    source_dir = sys.argv[2]
    database = sys.argv[3]
    alignment_dir = sys.argv[4]
    species_int = sys.argv[5]
    splitted_len = int(sys.argv[6])

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


        query_node = identifier + '_' + hit_dict["query"].replace('_+_','').replace('_-_','')
        ref_node = identifier + '_' + hit_dict["parent"]
        ancestor_node = identifier + '_' + hit_dict["ancestor"]

        bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                                    alignment_dir + identifier + '_pagan.fas',
                                    "fasta"))
        align_dict = OrderedDict()

        for seq in bio_align_dict:

            align_dict[identifier + '_' + seq.replace('_+_','').replace('_-_','')] = str(bio_align_dict[seq].seq)
            seq = identifier + '_' + seq.replace('_+_','').replace('_-_','')
            if 'mus_musculus' in seq:
                print(seq.split('_'))
            if species_int in seq and len(seq.split('_')) == splitted_len:

                species_int_node = seq

        source_start = hit_dict["tsm_source_start"] + 1
        source_end = hit_dict["tsm_source_end"]

        target_start = hit_dict["tsm_target_start"] + 1
        target_end = hit_dict["tsm_target_end"]

        for node in [query_node, ref_node, species_int_node]:
 
            filename = source_dir +  node + '.db'

            target_n = target_start
            
            if align_dict[node][target_start] == '-':
                target_n = target_start
                try:
                    while align_dict[node][target_n] == '-':
                        target_n += 1
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
                    '/')[-1]
            short_name = short_name_list + '_ss.ps'

            if node == query_node:
                tag = "qry"
            elif node == ref_node:
                tag = "ref"
            elif node == ancestor_node:
                tag = "anc"
            elif node == species_int_node:
                tag = species_int

            target_file = target_dir + identifier + '_' +\
                str(hit_dict["position"]) + "_" + tag + "_" +\
                short_name.rstrip('_ss.ps') + '.ps'
            shutil.copyfile(short_name.split(identifier+'_')[-1], target_file)

            if query_node == species_int_node:
                target_file = target_dir + identifier + '_' +\
                    str(hit_dict["position"]) + "_" + species_int + "_" +\
                    short_name.rstrip('_ss.ps') + '.ps'
                shutil.copyfile(short_name, target_file)


if __name__ == "__main__":
    main()
