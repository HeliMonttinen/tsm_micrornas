"""
This script collects information on those
gene region from a given .gff file, which overlap
with a microRNA. Writes an output file, in which
in the first column the gene identifier is given.
"""

import os
import sys
from collections import defaultdict

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def parse_gff_file(gff_file, mir_locs):
    """
    Parses a gff file and returns it as a dictionary.
    Chromosome and locations are keys. The value is
    a parsed line.

    Argement
    =========

    gff_file: A path pointing to a -gff file
    mir_locs: A dictionary containing microRNA
              locus.

    Returns
    ========

    gene_locations: A dictionary having chromosome and
                    location as a key. A line in a .gff file
                    is a value.

    """

    gene_locations = defaultdict(dict)

    with open(gff_file, 'r') as f:

        for line in f:

            line_splitted = line.rstrip().split()
            if len(line_splitted) >= 9:

                try:
                    chromosome = line_splitted[0]
                    start = int(line_splitted[3])
                    stop = int(line_splitted[4])
                    found = False

                    if chromosome in mir_locs:
                        for mir_start_stop in mir_locs[chromosome]:
                            mir_start = int(mir_start_stop.split('-')[0])
                            mir_end = int(mir_start_stop.split('-')[1])

                            if ((mir_end < start) or (mir_start > stop)):
                                continue
                            else:
                                found = True
                    if found is False:
                        continue
                    location = str(start) + '-' + str(stop)

                    if chromosome not in gene_locations:
                        gene_locations[chromosome] = dict()

                    if location not in gene_locations[chromosome]:
                        gene_locations[chromosome][location] = list()

                    gene_locations[chromosome][location].append(line.rstrip())

                except:
                    continue

    return gene_locations


def collect_mir_info(mir_file):
    """
    Collects microRNA identifiers, chromosome,
    and start_stop coordinates from a mir_file.

    Arguments
    =========

    mir_file: A .csv file containing microRNA
              identifiers and coordinates

    Returns
    =======

    miR_locs: A dictionary in which chromosome is
              the first key and the locus start_stop
              indexes as a second. The identier of a microRNA
              is a value.

    mir_identical:  dictionary in which chromosome is
              the first key and the locus start_stop
              indexes as a second. The value is
              a list of identifiers.
    """

    miR_locs = defaultdict(dict)
    mir_identical = defaultdict(dict)

    with open(mir_file, 'r') as f:

        for line in f:
            if 'Locus' not in line:

                line_splitted = line.split('\t')
                identifier = line_splitted[0]
                chromosome = line_splitted[2].split(':')[0]
                start_stop = line_splitted[2].split(':')[1]

                if chromosome not in miR_locs:
                    miR_locs[chromosome] = dict()

                if start_stop not in miR_locs[chromosome]:
                    miR_locs[chromosome][start_stop] = identifier
                else:
                    if start_stop not in mir_identical[chromosome]:
                        mir_identical[chromosome][start_stop] = list()
                    mir_identical[chromosome][start_stop].append(identifier)

    return miR_locs, mir_identical


def main():
    """
    Runs the code.

    Input
    =====
    mir_file: A path to a -csv file containing
              information on microRNAs of interest
    gff_file: A path to a gff file
    """

    mir_file = sys.argv[1]
    gff_file = sys.argv[2]

    mir_locs, mir_identical = collect_mir_info(mir_file)

    gene_locations = parse_gff_file(gff_file, mir_locs)

    for chrom in mir_locs:
        for locs in mir_locs[chrom]:
            identifier = mir_locs[chrom][locs]

            mir_ints = list()
            mir_ints.append(identifier)
            if identifier in mir_identical:
                for idx in mir_identical[identifier]:
                    mir_ints.append(idx)
            mir_ints = list(set(mir_ints))
            found = False
            miR_start = int((locs).split('-')[0])
            miR_stop = int((locs).split('-')[1])

            for gen_loc in gene_locations[chrom]:

                for item in  gene_locations[chrom][gen_loc]:
                    gen_loc_split = item.split()

                    if gen_loc_split[2] == "chromosome":
                        continue

                    try:
                        gene_start = int((gen_loc).split('-')[0])
                        gene_stop = int((gen_loc).split('-')[1])

                    except:
                        continue

                    if ((miR_stop < gene_start) or (miR_start > gene_stop)):
                        continue
                    else:

                        gen_start_stop = set(range(gene_start, gene_stop))

                        found = True

                        gene_length_one_third = int(len(gen_start_stop) / 3)

                        gene_begin = gene_start + gene_length_one_third
                        gene_mid = gene_begin + gene_length_one_third
                        start_overlap = False
                        mid_overlap = False
                        end_overlap = False

                        if miR_stop < gene_begin:
                            start_overlap = True
                        elif miR_stop < gene_mid:
                            mid_overlap = True
                        else:
                            end_overlap = True

                        for idsz in mir_ints:
                            if start_overlap is True:
                                str1 = idsz + '\t' + locs + '\t' +\
                                        '\t' + "start" + '\t' +\
                                        item
                                print(str1)

                            elif mid_overlap is True:

                                str1 = idsz + '\t' + locs + '\t' +\
                                        "mid" + '\t' +\
                                        item

                                print(str1)

                            elif end_overlap is True:
                                str1 = idsz + '\t' + locs + '\t' +\
                                        "end" + '\t' +\
                                        item
                                print(str1)

            if found is False:
                for idsz in mir_ints:
                    str1 = idsz + '\t' + locs + '\t' +\
                        '\t' + "None" + '\t' + "Functionally unknown region"
                    print(str1)


if __name__ == "__main__":
    main()
