"""
These tools are intended for support usage of
fpa.

Author: Heli Mönttinen (0000-0003-2461-0690)
"""

from collections import defaultdict


def fpa_parsing(fpa_file):
    """
    This tool is for parsing a fpa file.
    It goes through a given fpa file and returns a defaultdict
    dictionary as a result.

    Parameters
    ----------

    :fpa_file: A FPA file that is an output of the fpa
    python interface. Should contain a file tag for
    the sequence identifiers that are aligned.

    Returns
    -------

    fpa_dictionary: A default dict formatted dictionary.
                    A key is a tuple of two sequence identifiers
                    (in sorted order).
                    The result dict contains
                    starting and ending indexes and sequences  for
                    a tsm source and target regions.

    """

    fpa_dict = defaultdict(list)
    fpa_dict_solutions = defaultdict(list)

    with open(fpa_file, 'r') as f:
        solve_start = False
        for line in f:
            if "#Sequences" in line:

                if solve_start is True:

                    fpa_dict_solutions[sequence_ids[0],
                                       sequence_ids[1]].append(new_line)
                solve_start = True
                new_line = ""
                case = {
                    "ref": "",
                    "query": "",
                    "start": "",
                    "end": "",
                    "sw_start": "",
                    "sw_end": "",
                    "seq_qry": "",
                    "seq_ref": ""}

                sequence_ids = line.rstrip().split()[-1].split(",")
                case["query"] = sequence_ids[1]
                case["ref"] = sequence_ids[0]
                sequence_ids.sort()

            elif len(line.split(",")) > 2:
                indexes = line.rstrip().split(',')

                case["start"] = indexes[5]
                case["end"] = indexes[9]
                case["sw_start"] = indexes[7]
                case["sw_end"] = indexes[8]

            elif "qry" in line and "|" in line:

                case["seq_qry"] = line.split("|")[1]

                fpa_dict[sequence_ids[0], sequence_ids[1]].append(case)

                case["end"] = int(case["start"]) + \
                    len(case["seq_qry"].replace("-", ""))

            elif "ref" in line and "|" in line:

                case["seq_ref"] = line.split("|")[1]

            new_line += line

        fpa_dict_solutions[sequence_ids[0], sequence_ids[1]].append(new_line)

    return fpa_dict, fpa_dict_solutions
