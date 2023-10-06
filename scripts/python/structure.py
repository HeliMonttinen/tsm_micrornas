"""
Functions for structural analysis.

"""

import os
from subprocess import check_output
from RNA_alignments import indexes_in_alignment
from RNA_loops import _join_dot_file_lines


def sequences_for_single_structural_analysis(directory,
                                             alignment_file,
                                             data_dict,
                                             *args):
    """
    Creates a new directory for a structural analysis.
    Picks uniq identifiers from foa file and creates separate
    fasta files for them.

    :directory: A target directory
    :fpa: A target fpa
    :fasta: A fasta file including all the sequences of the cluster

    """

    identifiers = []
    for arg in args:
        identifiers.append(arg)

    if not os.path.exists(directory):
        os.mkdir(directory)

    for identifier in identifiers:
        new_dict = {}
        try:
            sequence = str((data_dict[identifier]).seq).replace('-', '')
        except:
            sequence = str(data_dict[identifier]).replace('-', '')
        new_dict[identifier] = sequence
        if directory not in identifier:
            file_for_data_dict(new_dict,
                            directory + alignment_file +\
                            identifier + '.fasta')
        else:
            file_for_data_dict(new_dict,
                            alignment_file +\
                            identifier + '.fasta')


def file_for_data_dict(data_dict, outputfile):
    """
    Makes a file for a data dictionary.
    Ensures that sequences are not aligned and
    removes '-' at first from the sequence.

    Arguments
    ----------
    :data_dict: A data dictionary for sequences
    :filename: A name for the output file

    Output
    ------
    A fasta file
    """

    with open(outputfile, 'w') as f:
        for identifier in data_dict:
            f.write('>' + identifier + '\n')
            f.write(data_dict[identifier] + '\n')


def analyse_structure(fasta_file,
                      outfile,
                      alignmentcommand):
    """
    runs alignment program for a set of sequences.
    Takes a set of sequences as a dictionary format.

    :in_file: fasta file of structures to be aligned
    :outfile: outfile of structures to be aligned
    :alignmentcommand: alignment command
    """

    p = check_output(
            ["{alignmentcommand}".format(
                alignmentcommand=alignmentcommand),
             "{fasta_file}".format(fasta_file=fasta_file)])

    with open(outfile, 'w') as f:
        f.write(str(p, 'utf-8'))


def identify_alignment_area(structure_file, align_dict, index1, index2):
    """
    Identifies from the indexes the area that overlaps with
    the alignment.

    :structure_file: The file that has the predicted structure
    :alignment_file: The file of the original alignment

    """

    structure = _join_dot_file_lines(structure_file)

    with open(structure_file, 'r') as f:
        for line in f:
            if '>' in line:
                structure_id = line.lstrip('>').rstrip()

    alignment_seq = align_dict[structure_id]

    aligned_struct = ""

    b = 0
    for a in range(len(alignment_seq)):
        if alignment_seq[a] == '-':

            aligned_struct = aligned_struct + '-'
        elif b < len(structure):
            aligned_struct = aligned_struct + structure[b]
            b += 1

    ts_region = aligned_struct[index1:index2]

    return ts_region, aligned_struct
