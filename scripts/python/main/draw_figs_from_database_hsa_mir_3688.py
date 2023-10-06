#!/usr/bin/env python
import os
import sys
from Bio import SeqIO
from ete3 import Tree
from collections import (defaultdict,
                         OrderedDict)
from pymongo import MongoClient

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """creates figures from ts cases that match the query."""

    from RNA_alignments import indexes_in_alignment
    from tree_parse import (fix_indexes_based_on_position,
                            subtree_for_visualization,
                            cut_align_seq_position)

    from structure import identify_alignment_area

    from visualize import (hit,
                           files,
                           params,
                           make_fig)

    output_dir = sys.argv[1]
    alignment_dir = sys.argv[2]
    tree_dir = sys.argv[3]
    structure_dir = sys.argv[4]
    microRNAs = sys.argv[5]
    micro_info = sys.argv[6]

    species_start_stop = {}
    attributes_dict = {}
    product_start_stop = defaultdict(dict)

    with open(micro_info, 'r') as f2:

        for line in f2:
            line_splitted = line.rstrip().split()
            if line_splitted[1] == "None":
                continue
            column_int = line_splitted[-1]
            mir_prods = column_int.split(';')
            identifier = line_splitted[0]
            mir_prod_dict = {}

            for prod in mir_prods:
                start = prod.split('(')[-1].split('-')[0]
                stop = prod.rstrip().rstrip(')').split('-')[-1]
                prod_name = prod.split('(')[0].lstrip('"').rstrip('"')
                
                if '-' in line_splitted[2]:
                    mir_prod_dict[prod_name] = {"start": int(start),
                                                "stop": int(stop)+2}
                else:
                    mir_prod_dict[prod_name] = {"start": int(start)-1,
                                                "stop": int(stop)+1}

            product_start_stop[identifier] = mir_prod_dict

    with open(microRNAs, 'r') as f:

        for line in f:
            line_splitted = line.split()
            
            locus = line_splitted[-1]
            if (locus.rstrip()).endswith('_'):
                continue
            species = '_'.join((locus.split('_'))[:2])
            start = int((locus.split(':')[-1]).split('-')[-2])
            end = int((locus.split(':')[-1]).split('-')[-1])

            ids = line_splitted[0]
            attributes_dict[locus] = ids
            species_start_stop[ids] = [start, end+1]
                
    used=set()


    hit_dict = {"identifier" : "ENSG00000264105",
                "miR_id" : "hsa-mir-3688",
                "position" : 0, 
                "human_chromosome" : "4",
                "human_location" : "4:159128752-159129043",
                "query" : "#5#",
                "parent" : "#9#",
                "sister" : "#8#",
                "ancestor" : "#29#",
                "query_child1" : "#2#",
                "query_child2" : "#4#",
                "tsm_source_start" : 67,
                "tsm_source_end" : 114,
                "tsm_target_start" : 133,
                "tsm_target_end" : 170,
                "quality" : True}

    identifier = hit_dict["identifier"]
    if identifier not in product_start_stop:

        mir_prod_dict["None"] = {"start": 10,
                                 "stop": 20}

        product_start_stop[identifier] = mir_prod_dict


    tree = Tree(tree_dir + identifier + '_pagan.anctree', format=1)

    bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                                alignment_dir + identifier + '_pagan.fas',
                                "fasta"))
    align_dict = OrderedDict()

    for seq in bio_align_dict:

        align_dict[seq] = str(bio_align_dict[seq].seq)

    position_0 = min(set([hit_dict["tsm_target_start"],
                          hit_dict["tsm_target_end"],
                          hit_dict["tsm_source_start"],
                          hit_dict["tsm_source_end"]]))

    position_1 = max(set([hit_dict["tsm_target_start"],
                          hit_dict["tsm_target_end"],
                          hit_dict["tsm_source_start"],
                          hit_dict["tsm_source_end"]]))

    orig_root = tree.get_tree_root()
    ancestor_node = orig_root.name
    ref_node = hit_dict["parent"]
    query_node = hit_dict["query"]
    parent = (tree&ref_node).up
    parent_name = parent.name
    qry_ts_start = hit_dict["tsm_target_start"]
    qry_ts_end = hit_dict["tsm_target_end"]
    ref_ts_start = hit_dict["tsm_source_start"]
    ref_ts_end = hit_dict["tsm_source_end"]
    
    used_loc = ref_node + ',' + query_node + ',' + str(qry_ts_start) + ',' + str(qry_ts_end)
    used.add(used_loc)

    minimum = min(set([qry_ts_start,
                       qry_ts_end,
                       ref_ts_start,
                       ref_ts_end]))

    maximum = max(set([qry_ts_start,
                       qry_ts_end,
                       ref_ts_start,
                       ref_ts_end]))

    query_node = hit_dict["query"]
    ref_node = hit_dict["parent"]
    sister_node = hit_dict["sister"]

    if "query_child1" in hit_dict:

        qry_child_structure_node1 = hit_dict["query_child1"]
        structure_file_qry_child1 = structure_dir + identifier + '_' +\
            hit_dict["query_child1"] + '.db'
    else:
        qry_child_structure_node1 = None
        structure_file_qry_child1 = None

    if "query_child2" in hit_dict:

        qry_child_structure_node2 = hit_dict["query_child2"]
        structure_file_qry_child2 = structure_dir + identifier + '_' +\
            hit_dict["query_child2"] + '.db'

    else:
        qry_child_structure_node2 = None
        structure_file_qry_child2 = None

    structure_file_ref = structure_dir + identifier + '_' +\
        hit_dict["parent"] + '.db'
    structure_file_qry = structure_dir + identifier + '_' +\
        hit_dict["query"] + '.db'
    structure_file_sister = structure_dir + identifier + '_' +\
        hit_dict["sister"] + '.db'
    structure_file_parent = structure_dir + identifier + '_' +\
        hit_dict["ancestor"] + '.db'

    mir_acc = defaultdict(dict)
    product_start_stop2 = defaultdict(dict)
    mir_acc["MIMAT0019223"] = {"start": 79,
                               "stop": 112}

    mir_acc["MIMAT0018116"] = {"start": 135,
                               "stop": 157}

    product_start_stop2[identifier] = mir_acc

    seq_int = "Homo_sapiens_chr4_159128758_159129042" 
    structure_file_int = structure_dir + identifier + '_' +\
        seq_int + '.db'

    index_start = 67
    index_end = 171
    fig_start = 0
    fig_end = len(align_dict[seq_int]) -2

    if len(align_dict[seq_int]) > index_end + 50:
        fig_end = fig_end
    else:
        fig_end = len(align_dict[seq_int])-2

    if index_start > 50:
        fig_start = 0
    else:
        fig_start = 0

    (minimum,
     maximum,
     cut_start,
     cut_end,
     cut_align_dict) = cut_align_seq_position(
     align_dict,
     fig_start,
     fig_start,
     fig_end,
     fig_end)

    ts_region_qry, aligned_struct_qry = identify_alignment_area(
            structure_file_qry,
            align_dict,
            cut_start,
            cut_end)

    ts_region_ref, aligned_struct_ref = identify_alignment_area(
            structure_file_ref,
            align_dict,
            cut_start,
            cut_end)

    if not (tree&query_node).is_leaf():

        showseq = [parent_name,
                   ref_node,
                   query_node,
                   sister_node,
                   qry_child_structure_node1,
                   qry_child_structure_node2]

        showstruct = [parent_name,
                      ref_node,
                      query_node,
                      sister_node,
                      seq_int,
                      qry_child_structure_node1,
                      qry_child_structure_node2]

        struct_files = {
                parent_name: structure_file_parent,
                ref_node: structure_file_ref,
                query_node: structure_file_qry,
                sister_node: structure_file_sister,
                seq_int: structure_file_int,
                qry_child_structure_node1: structure_file_qry_child1,
                qry_child_structure_node2: structure_file_qry_child2}

    else:

        showseq = [parent_name, ref_node, query_node, sister_node, seq_int]
        showstruct = [parent_name, ref_node, query_node, sister_node, seq_int]

        struct_files = {ref_node: structure_file_ref,
                        query_node: structure_file_qry,
                        sister_node: structure_file_sister,
                        parent_name: structure_file_parent,
                        seq_int: structure_file_int}

    if qry_ts_end > len(aligned_struct_ref):
        qry_ts_end = len(aligned_struct_ref)
        ref_ts_end = len(aligned_struct_ref)

    hit_class = hit(10,
                    300,
                    ref_node,
                    query_node,
                    qry_ts_start,
                    qry_ts_end,
                    ref_ts_start,
                    ref_ts_end,
                    aligned_struct_ref,
                    aligned_struct_qry,
                    showseq,
                    showstruct,
                    hit_dict["human_location"].split(':')[0],
                    hit_dict["human_location"].split(':')[1],
                    hit_dict["miR_id"],
                    hit_dict["identifier"],
                    67,
                    171,
                    product_start_stop2[identifier])

    par = params()

    par.lmar = 350

    files_class = files(
            align_dict,
            output_dir + identifier + '_' +
            str(hit_dict["position"]) + ".png",
            struct_files)

    make_fig(files_class,
             hit_class,
             par,
             tree)


if __name__ == "__main__":
    main()
