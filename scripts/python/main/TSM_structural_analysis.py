"""
This script goes through all identified
TSMs in gene regions that do not overlap
with a simple repeat. It predicts secondary structure
and classifies structural changes caused by a TSM.
In addition identifies, where in a human
lineage a TSM has happened. The output is written in
a .csv.
"""

from collections import (OrderedDict,
                         defaultdict)
import os
import sys
from Bio import SeqIO
from ete3 import Tree
from pymongo import MongoClient
import re

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def structure(alignment_dir,
              structure_dir,
              separate):
    """
    Predicts secondary structure for a sequence between
    minimum and maximum TSM event indexes.

    Arguments
    ==========

    alignment_dir: A directory where alignment files are saved.
    structure_dir: A directory where structural files are saved.
    separate: a dictionary where gene chunk and TSM number are keys
              and tree nodes as list values.
    """

    from structure import (sequences_for_single_structural_analysis,
                           analyse_structure)
    for identifier in separate:

        identifier_file = '_'.join(identifier.split('_')[:2])
        bio_align_dict = SeqIO.to_dict(
                SeqIO.parse(
                    alignment_dir + identifier_file + '_pagan.fas',
                    "fasta"))

        align_dict = OrderedDict()

        struct_list = []
        try:
            align_dict[separate[identifier][4]] =\
                str(bio_align_dict[separate[identifier][4]].seq)[
                    separate[identifier][0]-20:separate[identifier][3]+20]

        except:
            align_dict[separate[identifier][4]] =\
                    str(bio_align_dict[separate[identifier][4]].seq)[
                        separate[identifier][0]:separate[identifier][3]]
        try:
            align_dict[separate[identifier][5]] =\
                str(bio_align_dict[separate[identifier][5]].seq)[
                    separate[identifier][0]-20:separate[identifier][3]+20]
        except:
            align_dict[separate[identifier][5]] =\
                str(bio_align_dict[separate[identifier][5]].seq)[
                    separate[identifier][0]:separate[identifier][3]]

        if not os.path.exists(
                structure_dir + identifier + '_' +
                separate[identifier][4] + '.db'):
            struct_list.append(separate[identifier][4])
        if not os.path.exists(
                structure_dir + identifier + '_' +
                separate[identifier][5] + '.db'):
            struct_list.append(separate[identifier][5])

        sequences_for_single_structural_analysis(structure_dir,
                                                 identifier,
                                                 align_dict,
                                                 *struct_list)

        for struct in struct_list:
            try:

                analyse_structure(
                        structure_dir + identifier + '_' +
                        struct + ".fasta",
                        structure_dir + identifier + '_' +
                        struct + ".db",
                        "RNAfold")
            except:
                continue


def create_annotation_dicts(annotation_dir):
    """
    Create dictionaries from annotation files.
    Identifies annotation files based on the filename.

    Arguments
    =========

    annotation_dir: A directory where bed files for different
                    properties are located.

    Returns
    ========

    Nested dictionaries for different properties. A chromosome
    is the key of the upper layer dictionary, a starting index
    is the key for the inner dictionary.
    """

    CDS_dict = {}
    exon_dict = {}
    UTR3_dict = {}
    UTR5_dict = {}
    gene_dict = {}
    ncRNA_dict = {}

    for filename in os.listdir(annotation_dir):

        filen = os.path.join(annotation_dir, filename)

        with open(filen, 'r') as f:

            for line in f:

                line_splitted = line.rstrip().split('\t')

                chromosome = line_splitted[0]
                start = int(line_splitted[1])
                end = int(line_splitted[2])

                if chromosome not in CDS_dict:
                    CDS_dict[chromosome] = {}
                    exon_dict[chromosome] = {}
                    UTR3_dict[chromosome] = {}
                    UTR5_dict[chromosome] = {}
                    gene_dict[chromosome] = {}
                    ncRNA_dict[chromosome] = {}

                if "CDS" in filen:
                    CDS_dict[chromosome][start] = end
                elif "exon" in filen:
                    exon_dict[chromosome][start] = end
                elif "UTR3" in filen:
                    UTR3_dict[chromosome][start] = end
                elif "UTR5" in filen:
                    UTR5_dict[chromosome][start] = end
                elif "ncRNA" in filen:
                    ncRNA_dict[chromosome][start] = end
                elif "gene_merge" in filen:
                    gene_dict[chromosome][start] = end

    return CDS_dict, exon_dict, UTR3_dict, UTR5_dict, gene_dict, ncRNA_dict


def identify_overlap(database,
                     collection,
                     TSM_dict,
                     tree,
                     gene_coordinates,
                     ncRNA_coordinates,
                     CDS_coordinates,
                     exon_coordinates,
                     UTR3_coordinates,
                     UTR5_coordinates):
    """
    Identifies overlap between a TSM gene part or gene type.
    Identifies also, if TSM is an inplace inversion or not.

    Arguments
    ==========

    database: A mongo database where TSMs are saved.
    collection: A mongodatabase collection.
    TSM_dict: A dictionary about non-repeat TSMs
    tree: A TSM tree
    gene_coordinates: a dictionary about gene coordinates
    ncRNA_coordinates: a dictionary about ncRNA coordinates
    CDS_coordinates: a dictionary about CDS coordinates
    exon_coordinates: a dictionary about exon coordinates
    UTR3_coordinates: a dictionary about UTR3 coordinates
    UTR5_coordinates: a dictionary about UTR5 coordinates

    Returns
    ========
    inplace_inversion: A dictionary containing cases that are
                       classified as inplace inversions. They key is
                       TSM_chunk_idenfier + order number of a TSM.
                       The value is a list containing all TSM points
                       1-4 in a sorted order.
    separate: A dictionary containing cases that are
                       classified as inplace inversions. They key is
                       TSM_chunk_idenfier + order number of a TSM.
                       The value is a list containing all TSM points
                       1-4 in a sorted order.
    """

    client = MongoClient()
    db = client[database]
    collection = db[collection]

    inplace_inversion = defaultdict(list)
    separate = defaultdict(list)

    query = {"mismatches": {"$gt": 5},
             "TSM_length": {"$gt": 15},
             "quality": True}

    cursor = collection.find(
            query, no_cursor_timeout=True, batch_size=100)

    hit_dict_all = {}

    c = 0
    no_count = 0
    for document in cursor:

        hit_dict = {}
        c += 1
        for key, value in document.items():

            hit_dict[key] = value
        hit_dict_all[c] = hit_dict

    for a in hit_dict_all:

        hit_dict = hit_dict_all[a]

        source_start = hit_dict["tsm_source_start"]
        source_end = hit_dict["tsm_source_end"]
        target_start = hit_dict["tsm_target_start"]
        target_end = hit_dict["tsm_target_end"]

        human_location_start = int(hit_dict["human_location"].split('-')[0])

        TSM_identifier = hit_dict["human_chromosome"] + '\t' +\
            str(human_location_start + target_start) + '\t' +\
            str(human_location_start + target_end)

        if TSM_identifier not in TSM_dict:
            no_count += 1
            continue
        else:
            original_coordinates = TSM_dict[TSM_identifier]

        origin_start = int(original_coordinates.split(':')[-1].split('-')[0])
        origin_end = int(original_coordinates.split('-')[-1])

        (gene_type,
         gene_part) = identify_gene_part(
                 hit_dict["human_chromosome"],
                 origin_start,
                 origin_end,
                 gene_coordinates,
                 ncRNA_coordinates,
                 CDS_coordinates,
                 exon_coordinates,
                 UTR3_coordinates,
                 UTR5_coordinates)

        sort = sorted([source_start, source_end,
                       target_start, target_end])
        sort.append(hit_dict["query"])
        sort.append(hit_dict["parent"])
        sort.append(original_coordinates)
        sort.append(gene_type)
        sort.append(gene_part)
        common_ancestor = identify_location_in_the_tree(
                tree,
                hit_dict["taxonomic_group"])
        sort.append(common_ancestor)

        if (abs(source_start-target_start)) < 5 and\
                (abs(source_end-target_end)) < 5:

            inplace_inversion[
                hit_dict["identifier"] +
                '_' + str(hit_dict["position"])] =\
                        sort

            sort.append("in_place_inversion")
            separate[hit_dict["identifier"] +
                     '_' + str(hit_dict["position"])].extend(sort)

        else:

            sort.append("None")
            separate[hit_dict["identifier"] +
                     '_' + str(hit_dict["position"])].extend(sort)

    return inplace_inversion, separate


def identify_if_TSM_repeat(TSM_info):
    """
    This function parses a file containing
    information on TSM cases, which overlap with a simple
    repeat less than 50%. Creates a dictionary about coordiantes
    where coordinates in the sequence alignment are keys
    and the original coordiantes the value.

    Arguments
    ==========
    TSM_info: A path to a TSM_info .csv file

    Returns
    ========
    TSM_dict: A dictionary, where the TSM coordinates in the
    sequence alignment are keys and the original coordinates
    are values.
    """

    TSM_dict = {}

    with open(TSM_info, 'r') as f:

        for line in f:
            line = line.rstrip()
            line_split = line.split('\t')

            request = line_split[0] + '\t' + line_split[3] +\
                '\t' + line_split[4]
            original_coordinates = line_split[0] + ':' +\
                line_split[1] + '-' + line_split[2]

            TSM_dict[request] = original_coordinates
    return TSM_dict


def identify_gene_part(chromosome,
                       start,
                       end,
                       gene_dict,
                       ncRNA_dict,
                       CDS_dict,
                       exon_dict,
                       UTR3_dict,
                       UTR5_dict):
    """
    Go through the coordinates and identify the gene type and
    gene part where the TSM is located.

    Arguments
    ==========

    chromsome: chromosome (ucsc format)
    start: a starting index of the TSM insertion
    end: an ending index of the TSM insertion
    gene_dict: A path to a file containing gene coordinates
               (chr1 45435436    45437534)
    ncRNA_dict: A path to a file containing ncRNA gene coordinates
    CDS_dict: A path to a file containing CDS coordinates
    exon_dict: A path to a file containing exon coordinates
    3UTR_dict: A path to a file containing 3UTR coordinates
    5UTR_dict: A path to a file containing 5UTR coordinates


    Returns
    ========

    gene_type: A gene type (string)
    gene_part: A gene part (string)
    """

    found_in_gene = False
    found_in_gene_part = False
    gene_type = ""
    gene_part = ""

    for start_back in gene_dict[chromosome]:

        end_back = gene_dict[chromosome][start_back]

        if start_back <= start and end <= end_back:
            gene_type = "gene"
            found_in_gene = True

        elif start <= start_back and start_back < end:

            gene_type = "gene"
            found_in_gene = True

        elif start < end_back and end_back <= end:

            gene_type = "gene"
            found_in_gene = True

        elif start <= start_back and end_back <= end:
            gene_type = "gene"
            found_in_gene = True

    if found_in_gene is False:

        for start_back in ncRNA_dict[chromosome]:

            end_back = ncRNA_dict[chromosome][start_back]

            if start_back <= start and end <= end_back:
                gene_type = "ncRNA"
                found_in_gene = True

            elif start <= start_back and start_back < end:

                gene_type = "ncRNA"
                found_in_gene = True

            elif start < end_back and end_back <= end:

                gene_type = "ncRNA"
                found_in_gene = True

            elif start <= start_back and end_back <= end:
                gene_type = "ncRNA"
                found_in_gene = True

    if found_in_gene is True and found_in_gene_part is False:

        for start_back in CDS_dict[chromosome]:

            end_back = CDS_dict[chromosome][start_back]

            if start_back <= start and end <= end_back:
                gene_part = "CDS"
                found_in_gene_part = True

            elif start <= start_back and start_back < end:

                gene_part = "CDS"
                found_in_gene_part = True

            elif start < end_back and end_back <= end:

                gene_part = "CDS"
                found_in_gene_part = True

            elif start <= start_back and end_back <= end:
                gene_part = "CDS"
                found_in_gene_part = True

    if found_in_gene is True and found_in_gene_part is False:

        for start_back in exon_dict[chromosome]:

            end_back = exon_dict[chromosome][start_back]

            if start_back <= start and end <= end_back:
                gene_part = "exon"
                found_in_gene_part = True

            elif start <= start_back and start_back < end:

                gene_part = "exon"
                found_in_gene_part = True

            elif start < end_back and end_back <= end:

                gene_part = "exon"
                found_in_gene_part = True

            elif start <= start_back and end_back <= end:
                gene_part = "exon"
                found_in_gene_part = True

    if found_in_gene is True and found_in_gene_part is False:

        for start_back in UTR3_dict[chromosome]:

            end_back = UTR3_dict[chromosome][start_back]

            if start_back <= start and end <= end_back:
                gene_part = "UTR3"
                found_in_gene_part = True

            elif start <= start_back and start_back < end:

                gene_part = "UTR3"
                found_in_gene_part = True

            elif start < end_back and end_back <= end:

                gene_part = "UTR3"
                found_in_gene_part = True

            if start <= start_back and end_back <= end:
                gene_part = "UTR3"
                found_in_gene_part = True

    if found_in_gene is True and found_in_gene_part is False:

        for start_back in UTR5_dict[chromosome]:

            end_back = UTR5_dict[chromosome][start_back]

            if start_back <= start and end <= end_back:
                gene_part = "UTR5"
                found_in_gene_part = True

            elif start <= start_back and start_back < end:

                gene_part = "UTR5"
                found_in_gene_part = True

            elif start < end_back and end_back <= end:

                gene_part = "UTR5"
                found_in_gene_part = True

            elif start <= start_back and end_back <= end:
                gene_part = "UTR5"
                found_in_gene_part = True

    if found_in_gene is True and found_in_gene_part is False:

        gene_part = "intron"

    return (gene_type, gene_part)


def analyse_if_a_loop(separate,
                      structure_dir):
    """
    This function analyzes a structural change caused by
    a TSM. Hairpin structures are identified between the
    smallest and largest TSM coordinates.
    If the loop sequence is located in the middle third
    of the extracted sequence (between the smallest and
    largest TSM coordinates), the sequence is conseidered
    containing a hairpin. If both parent and child node
    pass this criteria, the hairpin is maintained.
    If only a parent node pass the criteria, a hairpin is
    destroyed. Id only a child node pass the criteria, a
    new hairpin is created.

    If stem length of a new hairpin is between 30 and
    55, the hairpin is considered as microRNA-like.

    Arguments
    ==========
    separate: An output of the function 'identify_overlap'.
              A dictionary where the key is a gene chunk identifier
              + a TSM number and the value is a list of
              TSM properties.
    structure_dir: A directory where secondary structure predioctions
                   of all TSM cases are saved.
    """

    new_hairpin = {}
    counts = {"new_hairpin": 0,
              "hairpin maintained": 0,
              "No strong structure": 0,
              "Hairpin destroyed": 0,
              "potential_microRNA": 0}
    information_dict = dict()

    for item in separate:

        identifier = '_'.join(item.split('_')[:2])
        position = item.rstrip().split('_')[-1]

        parent = False
        query = False
        potential_microRNA = False
        new_microRNA = False
        structural_change = "None"
        len_child = 0
        len_parent = 0
        for i in [separate[item][4], separate[item][5]]:

            if not os.path.exists(
                    structure_dir + item + '_' + i + '.db'):
                continue

            with open(structure_dir + item + '_' + i + '.db', 'r') as f:

                for line in f:

                    if '(' not in line and ')' not in line:

                        continue
                    else:
                        line = line.rstrip()
                        pattern = re.compile(r"(?<=\()[.-]+(?=\))")

                        longest_loop_start = 0
                        longest_loop_end = 0
                        longest_stem_start = 0
                        longest_stem_end = 0
                        longest_stem = 0

                        for loop in re.finditer(pattern, line.rstrip()):
                            indexes = loop.span()

                            loop_start = indexes[0]
                            loop_end = indexes[1]

                            stem_start = loop_start
                            stem_end = loop_end
                            longest_stem = 0

                            for a in reversed(line[:loop_start]):
                                if a == '(':
                                    stem_start -= 1
                                    stem_start_dot = 0
                                elif a == '.':

                                    stem_start_dot += 1
                                    stem_start -= 1

                                else:
                                    stem_start += stem_start_dot
                                    break

                            for a in line[loop_end:]:

                                if a == ')':
                                    stem_end += 1
                                    stem_end_dot = 0
                                elif a == '.':
                                    stem_end_dot += 1
                                    stem_end += 1
                                else:
                                    stem_end -= stem_end_dot
                                    break

                            if loop_start - stem_start > longest_stem:

                                longest_stem = loop_start - stem_start
                                longest_loop_start = loop_start
                                longest_loop_end = loop_end
                                longest_stem_start = stem_start
                                longest_stem_end = stem_end

                        if longest_stem == 0:
                            continue

                        if (len(line[
                            longest_stem_start:longest_loop_start])/(
                                longest_stem_end-longest_stem_start)) >\
                                0.33 and\
                                (len(line[
                                    longest_loop_end:longest_stem_end])/(
                                    longest_stem_end-longest_stem_start)) >\
                                0.33:

                            if i == separate[item][4]:
                                query = True
                                len_child = len(line[
                                    longest_stem_start:longest_loop_start])
                                if 30 <= len(line[
                                    longest_stem_start:longest_loop_start]) <=\
                                        55:
                                    potential_microRNA = True
                            else:
                                len_parent = len(line[
                                    longest_stem_start:longest_loop_start])
                                parent = True

        if (query is True and parent is False) or\
                (query is True and len_parent <= len_child - 5):
            new_hairpin[item] = "new_hairpin"
            counts["new_hairpin"] += 1
            structural_change = "new_hairpin"
            if potential_microRNA is True:
                counts["potential_microRNA"] += 1
                new_microRNA = True

        if query is True and parent is True and len_parent > len_child - 5:
            new_hairpin[item] = "hairpin maintained"
            counts["hairpin maintained"] += 1
            structural_change = "hairpin maintained"
        if query is False and parent is False:
            new_hairpin[item] = "No strong structure"
            counts["No strong structure"] += 1
            structural_change = "No strong structure"
        if query is False and parent is True:
            new_hairpin[item] = "Hairpin destroyed"
            counts["Hairpin destroyed"] += 1
            structural_change = "Hairpin destroyed"

        if identifier not in information_dict:
            information_dict[identifier] = defaultdict(list)
        separate[item].append(structural_change)
        separate[item].append(new_microRNA)
        information_dict[identifier][position].extend(separate[item])
    length = len(separate[item])
    return new_hairpin, counts, information_dict, length


def screen_secondary_TSMs(information_dict, length):
    """
    This function identifies, if a secondary TSMs have
    happened in the same locus with an older TSM event.

    Arguments
    ==========
    information_dict: An output of a function 'analyse_if_a_loop'.
                      It is a nested dictionary, in which
                      a gene chunk identifier is the key of the first
                      dictionary and TSM number the key of the
                      inner dictionary. The value is a list about TSM
                      properties.
    length: Length is a length of the TSM properties list.

    Returns
    ========
    Information_dict: The input dictionary, where added the information
                      whether the TSM is primary or secondary.
    """

    for item in information_dict:

        if len(information_dict[item].keys()) > 1:

            prev_coords = []
            for position in information_dict[item]:

                coord_pack = (int(information_dict[item][position][0]),
                              int(information_dict[item][position][3]),
                              int(information_dict[item][position][9]),
                              position)

                prev_coords.append(coord_pack)
            for coor in range(len(prev_coords)):
                for coor2 in range(len(prev_coords)):

                    if prev_coords[coor][2] < prev_coords[coor2][2]:

                        if prev_coords[coor][0] <=\
                                prev_coords[coor2][0] < prev_coords[coor][1]:

                            information_dict[item][
                                    prev_coords[coor][3]].append("secondary")
                        elif prev_coords[coor2][0] <=\
                                prev_coords[coor][0] < prev_coords[coor2][1]:

                            information_dict[item][
                                    prev_coords[coor][3]].append("secondary")

        for pos in information_dict[item]:

            if len(information_dict[item][pos]) == length:

                information_dict[item][pos].append("primary")

    return information_dict


def identify_location_in_the_tree(tree, species):
    """
    A fucntion to identify, where the TSM has happened
    in the phylogenetic tree.

    Arguments
    ==========

    tree: A phylogenetic tree
    species: A species list.

    Returns
    ========
    common_ancestor.name: A name of the node where the
                          TSM has happened.
    """

    species_present_nodes = []
    for spes in species:
        species_present_nodes.append(tree&spes)

    if len(species) == 1 and 'Homo' in species[0]:
        return 0

    common_ancestor = tree.get_common_ancestor(
                *species_present_nodes)

    return common_ancestor.name


def main():
    """
    Run the functions to analyze identified TSM
    cases.

    Arguments
    ==========

    TSM_info: A path to a file containing information on
              the TSMs not overlapping more than 50% with
              a simple repeat.
    database: a mongo database where TSM information is saved.
    collection: A mongo database collection
    structure_dir: A directory where structural files are saved.
    alignment_dir: A directory where alignment files are saved.
    tree_file: A path to a tree file.
    annotation_dir: A directory where gene coordinates are saved

    Output
    ======

    Collected information on TSM cases printed
    as tab-separated lines.
    """

    TSM_info = sys.argv[1]
    database = sys.argv[2]
    collection = sys.argv[3]
    structure_dir = sys.argv[4]
    alignment_dir = sys.argv[5]
    tree_file = sys.argv[6]
    annotation_dir = sys.argv[7]

    (CDS_dict,
     exon_dict,
     UTR3_dict,
     UTR5_dict,
     gene_dict,
     ncRNA_dict) = create_annotation_dicts(annotation_dir)

    TSM_dict = identify_if_TSM_repeat(TSM_info)

    tree = Tree(tree_file)
    tree_root = tree.get_tree_root()
    parent = tree&"Homo_sapiens"
    tree_node = tree&"Gorilla_gorilla_gorilla"
    tree_node.name = "Gorilla_gorilla"
    tree_node = tree&"Pongo_pygmaeus_abelii"
    tree_node.name = "Pongo_pygmaeus"

    count = 1
    while parent != tree_root:
        parent = (parent.up)
        parent.name = count
        count += 1

    (inplace_inversion,
     separate) = identify_overlap(
             database,
             collection,
             TSM_dict,
             tree,
             gene_dict,
             ncRNA_dict,
             CDS_dict,
             exon_dict,
             UTR3_dict,
             UTR5_dict)

    structure(alignment_dir,
              structure_dir,
              separate)

    new_hairpin, counts, information_dict, length = analyse_if_a_loop(
        separate,
        structure_dir)

    information_dict = screen_secondary_TSMs(information_dict, length)

    for item in information_dict:
        for position in information_dict[item]:
            string = item + '\t' + str(position) + '\t' +\
                    str(information_dict[item][position][0])
            for t in information_dict[item][position][1:]:
                string += '\t' + str(t).rstrip()
            print(string)


if __name__ == "__main__":
    main()
