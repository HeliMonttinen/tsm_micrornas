"""These scripts are intended for analysing RNA structure."""

from Bio.Seq import Seq
from itertools import islice
import common


def _join_dot_file_lines(dot_file):
    """
    Joins a dot parenthesis file (dot parenthesis part of a file)
    into one string.

    dot_file:  A filepath to dot parenthesis

    Returns
    dot_parenthesis_string: a dot parenthesis part of
    a file as one string.
    """

    dot_parenthesis_string = ""
    start_options = ['.', '(', '[', ')', ']', '{', '}', '<', '>']

    with open(dot_file, 'r') as data_file:
        for line in islice(data_file, 1, None):
            if line.startswith(tuple(start_options)):
                dot_parenthesis_string += line.rstrip().split()[0]

    return dot_parenthesis_string


def find_RNA_file(RNA_id, root_dir):
    """
    Returns a dot file, if a filename contains the given
    RNA identifier.

    RNA_id: RNA identifier
    root_dir: a root_dir to look_at

    Returns

    dot_file: a filepath to a dot file
    """

    for dot_file in common.return_file(root_dir):

        if RNA_id in dot_file:
            return dot_file

    return


def make_reverse_complement(seqs):
    """
    Takes a sequences (a list or string).
    and returns a reverse complement sequence as
    a string or a list dependening on the input.

    seqs: sequences as a string or a list of strings

    output: a reverse complement as a string or a list
    depending on the input

    if an input is not a string or a list it raises an error.
    """

    if isinstance(seqs, str):
        rna_seq = Seq(seqs)

        return rna_seq.reverse_complement()

    elif isinstance(seqs, list):

        reverse_rna_seqs = []

        for seq in seqs:
            rna_seq = Seq(seq, generic_rna)
            reverse_rna_seqs.append(rna_seq.reverse_complement())

        return reverse_rna_seqs

    else:
        raise TypeError(
                'Your Sequence cannot be recognized as a string or a list')
