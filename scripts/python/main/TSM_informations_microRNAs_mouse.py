"""
Collects information on the TSM cases, extracts
the source and target sequences of TSM cases and
runs quality control for them. All cases are
saved into a Mongo database.

See do_checks module for the quality control details.

As we are interested about TSMs creating human microRNA
loci. The query sequence between the TSM site minimum and maximum
indexes has to be less than length of the human sequence*1.1.
As otherwise, there have been notable deletions in human
microRNA loci after the TSM insertion.

"""

import os
import sys
from mongoengine import (connect,
                         disconnect)


dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    Goes through all identified TSM cases and performs quality control and
    finally saves cases with additional information to the
    to the mongodatabase. Cases that do not pass the quality
    criteria are also saved to the database, but their quality status is false.
    For the Mongo database structure, see TSM_tools for details.

    Arguments
    =========

    fpa_dir: A directory where the fpa output files are saved.
    database: A database where information on TSM cases is saved
    alignments_dir: A directory where the alignment files are located
    tree_dir: A direcotry where the tree files are located.
    mir_file: Information on the human microRNA loci
    fpa_solutions: A directory where FPA solutions are saved.
    species_int: The species, which miRNAs are studied. 

    Output
    =======
    A document from each TSM case in the Mongo database.

    """
    from collections import OrderedDict

    from Bio import SeqIO
    from ete3 import Tree

    from common import return_file
    from results_database import add_information

    from tree_parse import map_template_switch_in_same_region

    from TSM_tools import (sequence_information,
                           taxonomic_group,
                           miR_info)

    fpa_dir = sys.argv[1]
    database = sys.argv[2]
    alignments_dir = sys.argv[3]
    tree_dir = sys.argv[4]
    mir_file = sys.argv[5]
    fpa_solutions = sys.argv[6]
    species_int = sys.argv[7]

    connect(database, "default")

    mir_dict = miR_info(mir_file)

    used_locations = set()

    for filename in return_file(fpa_dir):

        identifier = filename.split('/')[-1]

        tree = Tree(tree_dir + identifier + '_pagan.anctree', format=1)

        bio_align_dict = SeqIO.to_dict(SeqIO.parse(
            alignments_dir + identifier + '_pagan.fas',
            "fasta"))

        align_dict = OrderedDict()

        for seq in bio_align_dict:
            align_dict[seq] = str(bio_align_dict[seq].seq)

        template_switch_clusters, solutions =\
            map_template_switch_in_same_region(
                    filename,
                    alignments_dir + identifier + '_pagan.fas',
                    do_not_merge=True)

        for position in template_switch_clusters:

            ref_node = template_switch_clusters[position][2][0]
            query_node = template_switch_clusters[position][2][1]


            qry_ts_start = int(template_switch_clusters[position][2][4])
            qry_ts_end = int(template_switch_clusters[position][2][5])
            ref_ts_start = int(template_switch_clusters[position][2][7])
            ref_ts_end = int(template_switch_clusters[position][2][6])

            used_location = ref_node + ',' + query_node + ',' +\
                str(qry_ts_start) + ',' + str(qry_ts_end)

            if used_location in used_locations:
                continue


            ref_tree = tree&template_switch_clusters[position][2][0]

            child_node1 = ref_tree.children[0].name
            child_node2 = ref_tree.children[1].name
            if child_node1 != template_switch_clusters[position][2][1]:

                sister_node = child_node1
            else:
                sister_node = child_node2

            align_ts_child1 = None
            align_ts_child2 = None
            qry_child1 = "None"
            qry_child2 = "None"

            if not (tree&query_node).is_leaf():

                qry_child1 = (tree&query_node).children[0].name
                qry_child2 = (tree&query_node).children[1].name
                align_ts_child1 = str(align_dict[qry_child1])\
                    [int(qry_ts_start):int(qry_ts_end)]
                align_ts_child2 = str(align_dict[qry_child2])\
                    [int(qry_ts_start):int(qry_ts_end)]

            align_ts_qry = str(align_dict[query_node])\
                [int(qry_ts_start):int(qry_ts_end)]
            align_ts_ref = str(align_dict[ref_node])\
                [int(qry_ts_start):int(qry_ts_end)]
            align_orig_ref = str(align_dict[ref_node])\
                [int(ref_ts_start):int(ref_ts_end)]
            align_orig_qry = str(align_dict[query_node])\
                [int(ref_ts_start):int(ref_ts_end)]

            minimum = min(set([qry_ts_start,
                              qry_ts_end,
                              ref_ts_end,
                              ref_ts_start]))
            maximum = max(set([qry_ts_start,
                               qry_ts_end,
                               ref_ts_end,
                               ref_ts_start]))

            quality = True

            (mismatches,
             length) = sequence_information(
                     align_ts_qry,
                     align_ts_ref)

            taxo_group = taxonomic_group(tree&query_node)
            
            for seq in align_dict:

                if species_int in seq and len(seq.split('_')) == 5:
                    species_int = seq
                    chromosome = seq.split('_')[2]
                    location = seq.split('_')[3] + '-'+ seq.split('_')[4]
                    species_int_seq =len(
                            (align_dict[seq][
                                minimum:maximum]).replace('-', ''))

            results_data = {
                    "identifier": identifier,
                    "miR_id": identifier,
                    "position": position,
                    "human_chromosome": chromosome,
                    "human_location": location,
                    "query": query_node,
                    "parent": ref_node,
                    "sister": sister_node,
                    "ancestor": "none",
                    "child1": qry_child1,
                    "child2": qry_child2,
                    "tsm_orig_start": ref_ts_start,
                    "tsm_orig_end": ref_ts_end,
                    "tsm_target_start": qry_ts_start,
                    "tsm_target_end": qry_ts_end,
                    "taxonomic_group": taxo_group,
                    "sequence_length": 15,
                    "TSM_length": length,
                    "mismatches": mismatches,
                    "quality": quality}

            used_locations.add(used_location)

            minimum = min(set(
                [ref_ts_start, ref_ts_end, qry_ts_start, qry_ts_end]))
            maximum = max(set(
                [ref_ts_start, ref_ts_end, qry_ts_start, qry_ts_end]))

            query_loc = len(
                    (align_dict[query_node][
                        minimum:maximum]).replace('-', ''))

            if '_'.join(species_int.split('_')[:2]) in taxo_group:
                
                add_information(results_data)
                with open(
                        fpa_solutions + identifier + '_' +
                        str(position) + '.txt', 'a+') as f3:
                    f3.write(solutions[position])

    disconnect("default")


if __name__ == "__main__":
    main()
