"""
Tools for TSM quality check and to
collect other information on TSMs.

"""

import os
import sys


dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def do_checks(align_dict,
              qry_ts_start,
              qry_ts_end,
              align_orig_ref,
              align_ts_ref,
              align_orig_sister,
              align_ts_sister,
              align_ts_qry,
              align_orig_qry,
              align_orig_anc,
              align_ts_anc,
              align_ts_anc_sis,
              align_orig_anc_sis,
              child_tar_seq1=None,
              child_tar_seq2=None):
    """
    Quality check for the sequence alignment quality in
    identified TSM region:

    Arguments
    ==========

    align_dict: A dictionary containing aligned sequences in str format.
    qry_ts_start: A starting index of the TSM target site in a sequence
                  alignment-
    qry_ts_end: The last index of the TSM target site in a sequence_alignment
    align_orig_ref: The sequence of a TSM source site in a parental node
    align_ts_ref: The sequence of a TSM target site in a parental node
    align_orig_sister: The sequence of a TSM source site in the sister
                       node of (a child node)
    align_ts_sister: The sequence of a TSM target siten in a sister node
    align_ts_qry: The sequence of a TSM target site in a child node
    align_orig_qry: The sequence of a TSM source site in a child node
    align_orig_anc: The sequence of a TSM source site in the ancestral
                    node of a parent
    align_ts_anc: The sequence of a TSM target site in the ancestral
                  node of a parent
    align_ts_anc_sis: The sequence of a TSM target site in the sister
                      of an ancestror node
    align_orig_anc_sis: The sequence of a TSM source site in the sister
                        of an ancestor node
    child_tar_seq1: The sequence of a TSM target siten in the child of
                    a child node (default:None)
    child_tar_seq2: The sequence of a TSM target siten in the child of
                    a child node (default:None)

    Returns
    =======
    """
    from tree_parse import (reference_sequence_quality,
                            compare_predicted_sequence)

    flag = False
    if reference_sequence_quality(
            align_orig_ref,
            align_ts_ref,
            align_orig_sister,
            align_ts_sister,
            align_orig_qry,
            ancestor_orig=align_orig_anc,
            ancestor_ts=align_ts_anc,
            anc_sis_ts=align_ts_anc_sis) is True:
        flag = True
    else:
        print("reference false")

    if compare_predicted_sequence(
            align_ts_anc,
            align_ts_anc_sis,
            relaxed=True,
            iupac=True) is False:
        flag = False

        print("anc_sis_ts_false")

    if compare_predicted_sequence(
            align_orig_anc,
            align_orig_anc_sis,
            relaxed=True,
            iupac=True) is False:
        flag = False

        print("anc_orig_false")
    # Quality check of the query sequence prediction, if
    # query is not a leaf
    if child_tar_seq1 is not None and flag is True:
        for child in [child_tar_seq1, child_tar_seq2]:
            if compare_predicted_sequence(
                    align_ts_qry,
                    set([child]),
                    iupac=True) is False:
                flag = False

            else:
                flag = True
                break

    if flag is False:
        return False

    return True


def sequence_information(align_ts_qry,
                         align_ts_ref):
    """
    Collects information on a sequence. E.g. mismatches,
    length, children.

    Argument
    =========

    align_ts_qry: The aligned sequence of a TSM target site
                  in a child node
    align_ts_ref; The aligned sequence of a TSM target site
                  in a parent node

    Returns
    =======

    mismathes: Number of mismatches between the given sequences
    length: The ungapped length of the TSM target site
            in the child sequence
    """
    from tree_parse import mismatches_between_sequences

    mismatches, length = mismatches_between_sequences(
            align_ts_qry,
            align_ts_ref)

    return mismatches, length


def taxonomic_group(subtree):
    """
    Returns taxonomic names in the given subtree

    Arguments
    ==========
    subtree: A subtree as an input

    Returns
    =======
    tax_names: A list of species included in a given subtree.
    """

    tax_names = set()
    for node in subtree.iter_leaves():

        tax_name_list = (node.name).split('_')
        tax_name = '_'.join(tax_name_list[:2])
        tax_names.add(tax_name)

    tax_names = list(tax_names)

    return tax_names


def miR_info(mir_file):
    """
    Reads a miRNA .csv file containing information on microRNAs
    and creates a dictionary.

    Arguments
    =========

    mir_file: A .csv file in contining colum Ensemble_identifier,
              _, mir name, location(chr:start-stop)

    Returns
    ========

    mir_dict: A dictionary containing microRNA information
              using the identifier in the first column
    """

    with open(mir_file) as f:

        mir_dict = {}

        for line in f:

            line_split = line.split()
            if len(line_split) == 4:
                mir_identifier = line_split[2]
            else:
                mir_identifier = "None"

            loc_last = int(line_split[-1].split(':')[1].split('-')[-1])
            loc_first = int(line_split[-1].split(':')[1].split('-')[0])
            info = {"mir_identifier": mir_identifier,
                    "chromosome": line_split[-1].split(':')[0],
                    "sequence": line_split[-1],
                    "length": int(loc_last - loc_first),
                    "loc_first": loc_first,
                    "loc_last": loc_last}

            mir_dict[line_split[0]] = info

    return mir_dict
