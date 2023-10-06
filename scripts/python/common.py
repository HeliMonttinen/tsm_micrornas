"""
A collection of general purpose scripts.

"""
import os


def return_file(root_dir):
    """
    Yields all the filepaths under a dir.
    root_dir: a directory to looked at

    Yields
    filepath: Path to the file.
    """

    for subdir, dirs, files in os.walk(root_dir):
        for filename in files:
            filepath = subdir + os.sep + filename

            yield filepath


def contains_any(string, chr_set):
    """
    Identifies, if a string contains any of the
    characters in a set.

    string: a string to look at.
    chr_set: a set of characters to search.

    Returns

    True if any of the charaters is found

    """

    return 1 in [char in string for char in chr_set]


def split_multiple_fasta_file(fasta, not_included=None,
                              include=None, unwanted=None):
    """
    Yields one fasta at a time. Not_included is an optional parameter,
    which includes identifiers that are not wanted to be
    returned.

    Requires:

    fasta: A full path to the multiple sequence fasta file

    not_included: A set of identifiers, which fasta
    sequences are not wanted to be returned.

    include: A set of identifiers, which should be included.
    unwanted: A set of chars that are not allowed to exists
              in the data set.

    """

    sequence = False
    fasta_sequence = ""
    identifier = ""

    with open(fasta, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>  File '):
                if sequence is True:
                    yield fasta_sequence

                identifier = line.split()[2].split('.')[0]
                sequence = True
                fasta_sequence = ""

            elif line.startswith('>'):
                if sequence is True:
                    yield fasta_sequence

                line = line.replace('> ', '>')
                identifier = line.split()[0].lstrip('>').rstrip().lstrip()
                sequence = True
                fasta_sequence = ""

            if not_included and sequence is True:
                if identifier in not_included:
                    sequence = False

            elif (include is not None) and sequence is True:

                if identifier not in include:
                    sequence = False

            if len(identifier) > 0 and len(fasta_sequence) == 0\
                    and sequence is True:
                fasta_sequence = '>' + identifier + '\n'

            if sequence is True and '>' not in line:
                if unwanted is None:
                    fasta_sequence += line.rstrip().strip()
                elif unwanted is not None:
                    if contains_any(line.rstrip().strip(), unwanted):
                        sequence = False
                        fasta_sequence = ""
                    else:
                        fasta_sequence += line.rstrip().strip()

            elif len(line) == 0:
                if sequence is True:
                    yield fasta_sequence
                sequence = False
                fasta_sequence = ""
                identifier = ""

        if sequence is True:
            yield fasta_sequence
