"""
Tools for parsing phylogenetic trees.

Author: Heli Mönttinen (ORCID: 0000-0003-2461-0690)

"""

from collections import (defaultdict,
                         OrderedDict)
from ete3 import Tree

from Bio import SeqIO

from fpa_tools import fpa_parsing
from RNA_alignments import (indexes_in_alignment,
                            parse_sequences_from_pairwise_mafft,
                            unusual_iupac_char)


def read_tree(tree_file):
    """
    Read a tree from a file,

    Requires

    tree_file: A full path to the tree file

    Returns

    a loaded tree

    """

    tree = Tree(tree_file, format=1)

    return tree


def go_through_branches(dictomy_tree):
    """
    This tool goes systematically through tree branches
    and samples systematically the tree by returning
    pairs of branches.

    Requires:

    tree: a tree from which polytomy is resolved

    Yields:
        Node name, identifier of sequence1 and identifier of sequence 2
    """
    node_number = 0

    for node in dictomy_tree.traverse("postorder"):

        if node.is_leaf() is not True:
            sample_seq1 = node.children[0].name
            sample_seq2 = node.children[1].name

            yield node.name, sample_seq1, sample_seq2

        node_number += 1


def rename_nodes_maf(tree_file, identifier, outfile,
                     alignment=False, align_outfile=None):
    """
    Unifies sequences headers and tree node names.

    As a result writes an alignment into a file.

    Arguments
    =========
    tree_dir: A directory, in which tree files are located
    identifir: Identifier which links alignment files and
               tree files.
    outfile: An outfile where sequences are written.
    """

    if alignment is not False:
        key_list = []
        fasta_dict = parse_sequences_from_pairwise_mafft(
                        alignment + identifier)

        for a in fasta_dict:
            key_list.append(a)

    tree = read_tree(tree_file)
    for node in tree.traverse("postorder"):

        if alignment is not False:

            if '_' in node.name:

                for a in fasta_dict:

                    if str(node.name) in a:
                        node.name = a

    tree.write(format=4,
               outfile=outfile)


def mismatches_between_sequences(seq1,
                                 seq2,
                                 strict_mode=None):
    """
    Returns the number of mismatches
    and the length of alignment between two sequences (does not
    take into account those positions which both have a gap).


    Arguments
    ==========

    Seq1: an aligned sequence 1
    Seq2: an aligned sequence 2
    strict_mode: if true, gaps are not calculated as mismatches
                 (default: None)

    Reuturns
    =========
    :mismatches: Number of  mismatches
    :length: Length of the TSM

    """

    mismatches = 0
    length = 0

    if strict_mode is True:
        for a, b in zip(seq1, seq2):

            if a != b and a != '-':

                mismatches += 1

            if a != '-':

                length += 1

        return mismatches, length

    for a, b in zip(seq1, seq2):

        if a != b:

            mismatches += 1

        if a != '-' or b != '-':
            length += 1

    return mismatches, length


def sequence_quality(seq1, seq2, relaxed=None, iupac=None):
    """
    The aim of this function is to confirm sequence
    quality of a predicted sequence.
    By default the mismatches between these
    two sequences should be less than 10 %.
    In the relaxed mode the mismatches has to be less than
    30%.  The script should return True, the sequences
    fulfil the requirement.

    :seq1: The first sequence to compare (string)
    :seq2: The second sequence to compare (string)
    :relaxed: If True, a relaxed mode is applied.
    :iupac: if True, the uncertain iupac chars in the same
            position of the alignment are interpreted as
            the same base between seq1 and seq2, if the the groups
            the possible bases overlap.

    Returns :
        Boolean
    """

    if iupac is not None:

        seq1, seq2 = unusual_iupac_char(seq1, seq2)

    mismatches, length = mismatches_between_sequences(seq1, seq2)

    if length <= 5:
        return True

    if relaxed is True:
        if (Decimal(mismatches)/Decimal(length)) <= 0.3:
            return True

    if (Decimal(mismatches)/Decimal(length)) <= 0.2:
        return True

    return False


def compare_predicted_sequence(seq,
                               children_seqs,
                               ancestor_seq=None,
                               relaxed=None,
                               iupac=None):
    """
    The aim is to confirm the quality of a predicted seq.
    Compares a sequence to the ancestral sequence
    and its possible child sequences. The sequence has to be always
    similar to the ancestral sequence, one of the child
    sequences.

    seq: The predicted sequence (aligned) which quality is wanted
         to check
    ancestor_seq: The ancestor sequence (aligned) of the seq
    children seqs: A set or a list of child sequences

    Returns:
    True if the sequence is of high quality
    """

    if ancestor_seq:
        if sequence_quality(
                seq, ancestor_seq, relaxed=True, iupac=iupac) is False:

            return False

    for child in children_seqs:
        if sequence_quality(seq, child, relaxed=relaxed, iupac=iupac) is False:
            return False

    return True


def reference_sequence_quality(reference_orig,
                               reference_ts,
                               sister_orig,
                               sister_ts,
                               qry_orig,
                               ancestor_orig=None,
                               ancestor_ts=None,
                               anc_sis_ts=None,
                               relaxed=None):

    """
    Ensures that ..
    1) Ts-region of reference sequence is similar to sister
       sequence and to the ancestral sequence

    2) Source ts-region has to be similar to the sister
       sequence or ts_sequence and ancestral sequence.

    Arguments
    ----------

    reference_orig: The origin sequence of the ts in the reference
    reference_ts: The template switch sequence in the reference
    sister_orig: The origin sequence in the sister sequence
    sister_ts: The template switch in the sister sequence
    ancestor_orig: The origin sequence in the ancestor of
    the reference ancestor_ts: The template switch sequence
    in the ancestor of reference sequence

    Output
    --------

    True if all the requirements are fulfilled. False, if some of them fails.

    """

    if compare_predicted_sequence(reference_ts,
                                  set([sister_ts]),
                                  ancestor_seq=ancestor_ts,
                                  relaxed=relaxed,
                                  iupac=True) is False:

        if anc_sis_ts:
            if compare_predicted_sequence(
                    reference_ts,
                    set([anc_sis_ts]),
                    ancestor_seq=ancestor_ts,
                    relaxed=True,
                    iupac=True) is False:

                return False
        else:
            return False

    if compare_predicted_sequence(reference_orig,
                                  set([qry_orig]),
                                  ancestor_seq=ancestor_orig,
                                  relaxed=relaxed,
                                  iupac=True) is False:

        if compare_predicted_sequence(reference_orig,
                                      set([sister_orig]),
                                      ancestor_seq=ancestor_orig,
                                      relaxed=relaxed,
                                      iupac=True) is False:

            return False

    return True


def location_of_template_switch(fpa_file, alignment_file):
    """
    The location of the template switch in the alignment.

    Parameters
    -------------
    :fpa_output: The output file of fpa
    :alignment: Alignment as a dictionary format.

    Output
    --------
    :indexes: of the template switch in the alignment

    """

    fpa_dictionary, fpa_dict_solution = fpa_parsing(fpa_file)
    align_dict = SeqIO.to_dict(SeqIO.parse(alignment_file, "fasta"))

    for hit in fpa_dictionary:
        number = 0

        for case in fpa_dictionary[hit]:

            if int(case["sw_end"]) > 0:

                template_switch_indexes = [int(case["start"]),
                                           int(case["start"]) +
                                           (int(case["sw_start"]) -
                                            int(case["sw_end"])+1),
                                           int(case["sw_start"]),
                                           int(case["sw_end"])-1]
            else:
                template_switch_indexes = [int(case["start"]),
                                           int(case["start"]) +
                                           (int(case["sw_start"]) -
                                               int(case["sw_end"])+1),
                                           int(case["sw_start"]),
                                           int(case["sw_end"])]

            alignment_indexes2 = indexes_in_alignment(
                template_switch_indexes[0:2], align_dict[case["query"]].seq)

            alignment_indexes = indexes_in_alignment(
                    template_switch_indexes[2:], align_dict[case["ref"]].seq)

            solution = fpa_dict_solution[hit][number]

            number += 1
            if len(alignment_indexes2) > 1 and len(alignment_indexes) > 1:
                yield (case["ref"],
                       case["query"],
                       case["seq_ref"],
                       case["seq_qry"],
                       alignment_indexes2[0],
                       alignment_indexes2[1],
                       alignment_indexes[0],
                       alignment_indexes[1],
                       solution)


def map_template_switch_in_same_region(fpa_file,
                                       alignment_file,
                                       do_not_merge=None):
    """
    Identifies if there are several template switch in the same
    region of the alignment. Returns a regions dictionary with region_numbers
    in which there are tuples of reference, query, sequence,
    ind1, ind2, ind3 and ind4.

    Requires
    --------

    :fpa_file: fpa_file output with sequence pair information
    :alignment_file: alignment file

    Returns
    -------

    :regions: Regions defaultdict list. Each key corresponds
              to one template switch region. A list contains
              information on the node pairs that has a template
              switch in the same region. The first index of the
              list tells the starting index of the region, and
              second the ending index of the region.

    """

    regions = defaultdict(list)
    solutions = {}

    template_switch_region = 0

    for ref, query, seq_qry, seq_ref, start,\
            end, sw_start, sw_end, solution in\
            location_of_template_switch(fpa_file, alignment_file):

        template_switch = set(list(range(start, end+1)))

        found = False

        if do_not_merge is True:
            minimum = int(start)
            maximum = int(end)

            regions[template_switch_region].append(minimum)
            regions[template_switch_region].append(maximum)
            regions[template_switch_region].append(
                (ref, query, seq_qry, seq_ref, start, end, sw_start, sw_end))

            solutions[template_switch_region] = solution
            template_switch_region += 1

        if do_not_merge is None:
            for region in regions:
                region_range = list(
                        range(regions[region][0], regions[region][1]+1))

                if len(template_switch.intersection(region_range)) > 0:

                    minimum = min(set(region_range).union(template_switch))
                    maximum = max(set(region_range).union(template_switch))

                    regions[region][0] = minimum
                    regions[region][1] = maximum

                    regions[region].append(
                            (ref,
                             query,
                             seq_qry,
                             seq_ref,
                             start,
                             end,
                             sw_start,
                             sw_end))

                    found = True
                    break

            if found is False:
                regions[template_switch_region].append(start)
                regions[template_switch_region].append(end)
                regions[template_switch_region].append(
                    (ref,
                     query,
                     seq_qry,
                     seq_ref,
                     start,
                     end,
                     sw_start,
                     sw_end))

                solutions[template_switch_region] = solution

                template_switch_region += 1

    return regions, solutions


def fix_indexes_based_on_position(position_list,
                                  index,
                                  is_list=False):
    """
    Takes a position list (designed especially for the subtree
    visualization) and substracts from given indexes the number of
    those positions that are smaller than the given index.

    :position_list: list of indexes that are removed form the original
                    alignment
    :index_list: Indexes that should be adjusted to the altered
                 alignment length.
    """

    smallerThan = lambda x, y: [i for i in x if i < y]

    if is_list is True:
        fixed_index = []
        for index1 in index:

            smaller = smallerThan(position_list, index1)

            new_index = index1 - len(smaller)

            fixed_index.append(new_index)
        return fixed_index

    else:
        smaller = smallerThan(position_list, index)
        fixed_index = index - len(smaller)

        return fixed_index


def subtree_for_visualization(tree,
                              ancestor_node,
                              alignment_dict):
    """
    This script is designed for large trees that are impossible
    to visualize if the wholse tree is used.
    :subtree: a sub tree of interest
    :alignment_dict: information of alignment dictonary
    :tree_leaves: optional parameter that identifies, how many leaves
                  a tree has to have that this script is run. default: 200

    Returns:
        new_align_dict: A dictionary of aligned sequences from which
                        all the positions
                        removed that are gaps in all the aligned sequences.
    """

    new_align_dict = OrderedDict()

    subtree = tree&ancestor_node
    names = set()
    for node in subtree.traverse():
        names.add(node.name)

    for node in alignment_dict:
        if node in names:
            new_align_dict[node] = alignment_dict[node]
            last = node

    positions = []
    for position in range(len(new_align_dict[last])):

        gap = True

        for seq in new_align_dict:

            if new_align_dict[seq][position] != '-':
                gap = False
                continue

        if gap is True:
            positions.append(position)

    positions.reverse()

    for pos in positions:

        for seq in new_align_dict:

            new_align_dict[seq] =\
                    new_align_dict[seq][0: pos:] +\
                    new_align_dict[seq][pos + 1::]

    return new_align_dict, positions


def cut_align_seq_position(align_dict,
                           qry_ts_start,
                           qry_ts_end,
                           ref_ts_start,
                           ref_ts_end):
    """
    Cut sequences to match the sequence region of interest.
    Leaves -+ 20 residues around the region of interest.
    Returns a dictionary containing the cut sequences.

    Arguments:
    ------------
    :align_dict: Alignment file in dictionary format
    :qry_ts_start: The starting index of the TSM target region
    :qry_ts_end: The ending index of the TSM target source
    :ref_ts_start: The starting index of the TSM source region
    :ref_ts_end: The encing index of the TSM source region

    Returns:
    --------

    :minimum: The smallest index within TSM source
              and target regions
    :maximum: The largest index within the TSM source and
              target regions
    :cut_start: The first index of the cut sequence
    :cut_end: The last index of the cut seqeucne
    :cut_align_dict: A dictionary containing the cut sequences
    """

    minimum = 6000
    maximum = 0

    num_set = set([qry_ts_start, qry_ts_end, ref_ts_start, ref_ts_end])

    if min(num_set) < minimum:
        minimum = min(num_set)

    if max(num_set) > maximum:
        maximum = max(num_set)

    cut_align_dict = {}
    cut_start = 0
    cut_end = 0

    for align in align_dict:

        new_name = align

        if minimum > 20 and (len(align_dict[align]) - maximum) > 20:
            cut_align_dict[new_name] =\
                    str(align_dict[align][minimum - 20:maximum + 20])

            if cut_end == 0:
                cut_start = minimum
                cut_end = maximum

        elif minimum <= 20 and (len(align_dict[align]) - maximum) <= 20:
            cut_align_dict[new_name] =\
                    str(align_dict[align][0:-1])

            if cut_end == 0:
                cut_start = 0
                cut_end = len(align_dict[align]) - 1

        elif minimum <= 20:
            cut_align_dict[new_name] =\
                str(align_dict[align][0:maximum + 20])

            if cut_end == 0:
                cut_start = 0
                cut_end = maximum

        elif (len(align_dict[align]) - maximum) <= 20:
            cut_align_dict[new_name] =\
                    str(align_dict[align][minimum-20:-1])

            if cut_end == 0:
                cut_start = minimum - 20
                cut_end = len(align_dict[align]) - 1

    return minimum, maximum, cut_start, cut_end, cut_align_dict
