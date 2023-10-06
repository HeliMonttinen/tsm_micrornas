"""
This module identifies transposable elements
nearby a microRNA and creates a tab-separated file
from them. The transposable elements that potentially
couöd provide a promoter for microRNA are marked with
overlap_xxxx or near_L1 (includes also SINEs). The
transposable elements, which cannot provide a promoter
or are more than 1000 bases away are marked with 'near'.
"""

from collections import defaultdict
import sys


def collect_gene_coordinates(filename):
    """
    Parses coordinates from a microRNA file.

    Arguments
    ========
    filename: A file about microRNA locations

    Returns
    ========
    coordinates: A nested dictionary in which a chromosome
                 is a key of the first dictionary
                 and start and stop indexes the keys of the
                 second one. The Ensembl identifier is
                 a value.
    """

    coordinates = defaultdict(dict)

    with open(filename, 'r') as f:
        for line in f:

            line_splitted = line.split()
            chromosome = 'chr' + line_splitted[1].split(':')[0]
            start = int(line_splitted[1].split(':')[1].split('-')[0])
            end = int(line_splitted[1].rstrip().split(':')[1].split('-')[1])
            coordinates[chromosome][(start, end)] = line_splitted[0].rstrip()
    
    return coordinates

def identify_repeats(filename, gene_coordinates,
                     transposons_dirs, strand_sense):
    """
    Identify the repeats (transposable elements and other repeats)
    within 300,000 bases from the microRNA. Gives also
    an estimate, if the transposable element or repeat overlaps
    with a microRNA or and if it can provide a promoter.

    If a repeat element is SINE, LINE or LTR and their status
    is not 'near', they potentially could provide a promoter
    for a microRNA locus.

    Arguments
    ==========
    filename: A path to a Repetamasker output file.
    gene_coordinates: microRNA coordinates as dictionry.
                      an output from a function gene coordinates.
    transposons_dirs: A directory where an output file is written.
    strand_sense: A dictionary, where an Ensembl identifier is a key
                  and a strandness is a value.

    Returns
    =======
    a .csv file written into a 'results_2' file into a
    transposons_dir directory.

    """

    with open(filename, 'r') as f1:

        for line in f1:
            if len(line.rstrip()) == 0:
                continue
            if 'class' in line or 'position' in line in line\
                    or 'matching' in line:
                continue

            line_splitted = line.split()
            chrom = line_splitted[4].rstrip()

            start = int(line_splitted[5].rstrip())
            end = int(line_splitted[6].rstrip())
            sense = line_splitted[8].rstrip()
            if sense == 'C':
                sense = '-'
            transposon = line_splitted[10].rstrip()


            for gene in gene_coordinates[chrom]:

                microRNA_sense = strand_sense[
                        gene_coordinates[chrom][gene]]

                if gene[0] <= start <= gene[1] <= end and\
                        '-' in microRNA_sense and\
                        sense == '-':

                    with open(transposons_dirs + 'transposons_microRNA.txt', 'a+') as f2:
                        
                        f2.write(gene_coordinates[chrom][gene] + '\t' +
                                 transposon + '\t' + 'overlap_partial_prom' +
                                 '\t' + line_splitted[8] + '\t' +
                                 chrom.lstrip('chr') + '\t' + str(start) +
                                 '\t' + str(end) + '\n')

                elif gene[0] <= start <= gene[1] <= end:

                    with open(transposons_dirs + 'transposons_microRNA.txt', 'a+') as f2:

                        f2.write(gene_coordinates[chrom][gene] + '\t' +
                                 transposon + '\t' + 'overlap_partial' +
                                 '\t' + line_splitted[8] + '\t' +
                                 chrom.lstrip('chr') + '\t' + str(start) +
                                 '\t' + str(end) + '\n')

                elif start <= gene[0] <= end <= gene[1] and\
                        '+' in microRNA_sense and\
                         sense == '+':

                    with open(transposons_dirs + 'transposons_microRNA.txt', 'a+') as f2:

                        f2.write(gene_coordinates[chrom][gene] + '\t' +
                                 transposon + '\t' + 'overlap_partial_prom' +
                                 '\t' + line_splitted[8] + '\t' +
                                 chrom.lstrip('chr') + '\t' + str(start) +
                                 '\t' + str(end) + '\n')

                elif start <= gene[0] <= end <= gene[1]:

                    with open(transposons_dirs + 'transposons_microRNA.txt', 'a+') as f2:

                        f2.write(gene_coordinates[chrom][gene] + '\t' +
                                 transposon + '\t' + 'overlap_partial' +
                                 '\t' + line_splitted[8] + '\t' +
                                 chrom.lstrip('chr') + '\t' + str(start) +
                                 '\t' + str(end) + '\n')


                elif start <= gene[0] and gene[1] <= end and\
                        sense in microRNA_sense:

                    with open(transposons_dirs + 'transposons_microRNA.txt', 'a+') as f2:

                        f2.write(gene_coordinates[chrom][gene] + '\t' +
                                 transposon + '\t' + 'overlap_fully_covered' +
                                 '\t' + line_splitted[8] + '\t' +
                                 chrom.lstrip('chr') + '\t' + str(start) +
                                 '\t' + str(end) + '\n')

                elif ((end <= gene[0] <= end + 1000) and
                        (('+' in microRNA_sense and sense == '+') or
                         ('L1' in transposon and '+' in microRNA_sense))) or\
                        ((start-1000 <= gene[1] <= start) and
                         (('-' in microRNA_sense and sense == '-') or
                          ('L1' in transposon and '-' in microRNA_sense))):

                    with open(transposons_dirs + 'transposons_microRNA.txt', 'a+') as f2:

                        f2.write(gene_coordinates[chrom][gene] + '\t' +
                                 transposon + '\t' + 'near_prom' + '\t' +
                                 line_splitted[8] + '\t' +
                                 chrom.lstrip('chr') + '\t' + str(start)
                                 + '\t' + str(end) + '\n')

                elif (end <= gene[0] <= end + 1000) or\
                        (start-1000 < gene[1] <= start):
                    with open(transposons_dirs + 'transposons_microRNA.txt', 'a+') as f2:
                        f2.write(gene_coordinates[chrom][gene] + '\t'
                                 + transposon + '\t' + 'within_1000' + '\t' +
                                 line_splitted[8] + '\t' +
                                 chrom.lstrip('chr') + '\t' + str(start) +
                                 '\t' + str(end) + '\n')

                elif (end <= gene[0] <= end + 3000000) or\
                        (start-3000000 < gene[1] <= start):
                    with open(transposons_dirs + 'transposons_microRNA.txt', 'a+') as f2:
                        f2.write(gene_coordinates[chrom][gene] + '\t'
                                 + transposon + '\t' + 'near' + '\t' +
                                 line_splitted[8] + '\t' +
                                 chrom.lstrip('chr') + '\t' + str(start) +
                                 '\t' + str(end) + '\n')


def main():
    """
    Run the functions for finding transposable elements
    nearby the microRNAs.

    Arguments
    ===========
    gene_coordinates: A path to a file containing microRNA gene coordinates
    transposons_file: A path to a Repeatmasker outputfile
    transpospns_results_dir: A path to a directory where
    'results_2' output file is written-
    TSM_strands: A path to a file containing information on
    the TSM containing microRNAs.
    """
    gene_coordinates = sys.argv[1]
    transposons_file = sys.argv[2]
    transposons_results_dir = sys.argv[3]
    TSM_strands = sys.argv[4]

    strand_dict = {}
    with open(TSM_strands, 'r') as f:

        for line in f:

            line_splitted = line.split('\t')
            strand_dict[line_splitted[0]] = line_splitted[2].rstrip()

    gene_dict = collect_gene_coordinates(gene_coordinates)

    identify_repeats(transposons_file,
                     gene_dict,
                     transposons_results_dir,
                     strand_dict)


if __name__ == "__main__":
    main()
