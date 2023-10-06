"""
Gets gene coordinates from a .gff file.
Merge the overlapping gene regions and divide the regions
into ~1000 bases long chunks. Adds 25 bases at both
ends so that each chunk overlaps 50 bases with adjacent
chunks.
"""

import sys


def get_coordinates(line):
    """
    Parses coordinates from a gff3line.

    Argument
    ========
    line: A gene line from a gff3 file.

    Returns
    =======

    parent: A Gene identifier
    chrom: A chromose
    Start: A starting coordinate of a gene
    stop: The last coordinate of a gene

    """

    line_splitted = line.split()
    chrom = line_splitted[0]
    start = line_splitted[3]
    stop = line_splitted[4]
    if "ID=gene:" in line:
        parent = line_splitted[8].split("ID=gene:")[-1].split(";")[0]
    else:
        parent = None
    return parent, chrom, start, stop


def main():
    """
    Identify the gene chunks and write the gene coordinate output file.


    Arguments
    ==========
    gff_file: A human hg38 .gff3 file
    exon_output: Where the gene chunk coordinates are written.
    """

    gff_file = sys.argv[1]
    gene_output = sys.argv[2]

    gene_dict = {}

    with open(gff_file, "r") as f:

        for line in f:

            if ("###" not in line) and len(line.split()) < 3:
                continue

            if (
                len(line.split()) > 7
                and "gene" in line.split()[2]
                and "pseudo" not in line.split()[2]
                and "mirbase" not in line
                and "segment" not in line
            ):

                (gen_parent,
                 gen_chrom,
                 gen_start,
                 gen_stop) = get_coordinates(line)
                if gen_parent is None:
                    continue

                gen_covered = False
                if gen_chrom in gene_dict:
                    for items in gene_dict[gen_chrom]:
                        if (items[0] <= int(gen_start)) and\
                                (int(gen_stop) <= items[1]):
                            gen_covered = True

                if gen_covered is True:
                    continue

                if gen_chrom not in gene_dict:
                    gene_dict[gen_chrom] = dict()

                covered = True

                with open(gene_output, "a+") as f1:
                    if int(gen_stop) - int(gen_start) < 2000:

                        f1.write(
                            gen_parent
                            + "_"
                            + str(1)
                            + "\t"
                            + gen_chrom
                            + ":"
                            + gen_start
                            + "-"
                            + gen_stop
                            + "\n"
                        )

                    else:
                        prev = int(gen_start)
                        number = 1
                        for a in range(int(gen_start) + 1000,
                                       int(gen_stop), 1000):
                            covered = False
                            if a + 1000 < int(gen_stop):
                                st = prev - 25
                                en = a + 25
                                for b in gene_dict[gen_chrom]:
                                    if (b[0] <= st) and (en <= b[1]):
                                        covered = True
                                        prev = a
                                        break
                                if covered is False:

                                    f1.write(
                                        gen_parent
                                        + "_"
                                        + str(number)
                                        + "\t"
                                        + gen_chrom
                                        + ":"
                                        + str(st)
                                        + "-"
                                        + str(en)
                                        + "\n"
                                    )
                                    number += 1
                                    prev = a
                            else:
                                st = prev - 25

                                for b in gene_dict[gen_chrom]:
                                    if (b[0] <= st) and\
                                            (int(gen_stop) <= b[1]):
                                        covered = True
                                        break
                                if covered is False:
                                    f1.write(
                                        gen_parent
                                        + "_"
                                        + str(number)
                                        + "\t"
                                        + gen_chrom
                                        + ":"
                                        + str(st)
                                        + "-"
                                        + gen_stop
                                        + "\n"
                                    )

                gene_dict[gen_chrom][(int(gen_start), int(gen_stop))] = ""


if __name__ == "__main__":
    main()
