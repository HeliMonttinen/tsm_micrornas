"""
Go through the mutation cluster coordinates
and TSM coordinates, identify their overlap and
return the coordinates of these overlaps.
"""

import sys


def read_variant_file(variant_file):
    """
    Go through the human variant cluster
    coordinates. and return them as dict.

    Argument
    =========

    variant_file: A path to human variant file
                  .csv file (chr    1   10)

    Reutrns
    ========

    variant_dict: A nested dictionary where a chromosome
                  is a key of the outer dictionary and
                  a starting index a key of the inner
                  dictionary. The last index is a value.
    """

    variant_dict = {}
    with open(variant_file, 'r') as f:
        for line in f:

            line_splitted = line.rstrip().split('\t')

            chromosome = line_splitted[0]
            if 'chr' not in line_splitted[0]:
                chromosome = 'chr' + line_splitted[0]

            if chromosome not in variant_dict:
                variant_dict[chromosome] = {}

            variant_dict[chromosome][int(line_splitted[1])] =\
                int(line_splitted[2])

    return variant_dict


def identify_variant_regions_in_TSM(TSM_file,
                                    variant_dict):
    """
    Identify overlaps between TSM regions and
    mutation clusters. Yields coordinates of
    the overlapping regions.

    Arguments
    ==========
    TSM file: A file containing TSM coordinates
    variant: An output of the function 'read_variant_file'

    Yields
    =======
    coordinates of an overlapping region.

    """

    with open(TSM_file, 'r') as f:

        for line in f:

            line_splitted = line.rstrip().split('\t')

            chromosome = 'chr' + line_splitted[0]
            start = int(line_splitted[1])
            end = int(line_splitted[2])

            found = False
            for coord in variant_dict[chromosome]:

                if start <= coord < variant_dict[chromosome][coord] <= end:
                    found = True

                    yield chromosome + '\t' + str(coord) + '\t' +\
                        str(variant_dict[chromosome][coord])

                elif coord <= start < end <= variant_dict[chromosome][coord]:
                    found = True
                    yield chromosome + '\t' + str(start) + '\t' + str(end)

                elif start <= coord <= end <= variant_dict[chromosome][coord]:
                    found = True
                    yield chromosome + '\t' + str(coord) + '\t' + str(end)

                elif coord <= start <= variant_dict[chromosome][coord] <= end:
                    found = True
                    yield chromosome + '\t' + str(start) + '\t' +\
                        str(variant_dict[chromosome][coord])

                if found is True:
                    break

            if found is False:
                print(chromosome, start, end)


def main():
    """
    Run the functions.

    Arguments
    ==========

    TSM_file: A file where TSM coordinates are saved as .csv.
    variant_file: A file where human variant cluster coordinates
                  are saved as a .csv format.
    output_file: A file where the overlap coordinates are saved.
    """

    TSM_file = sys.argv[1]
    variant_file = sys.argv[2]
    output_file = sys.argv[3]

    variant_dict = read_variant_file(variant_file)

    for coordinate in identify_variant_regions_in_TSM(
            TSM_file,
            variant_dict):

        with open(output_file, 'a+') as f5:

            f5.write(coordinate + '\n')


if __name__ == "__main__":
    main()
