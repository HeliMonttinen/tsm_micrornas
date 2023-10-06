"""
This script goes through the sample file
and randomly created background files and statistical
analysis.
"""

import os
import sys

from decimal import Decimal


def create_annotation_dicts(annotation_dir):
    """
    Create dictionaries from annotation files.
    Identifies annotation files based on the filename.

    Arguments
    =========

    annotation_dir: A directory where bed files for different
                    properties are located.

    Returns
    ========

    Nested dictionaries for different properties. A chromosome
    is the key of the upper layer dictionary, a starting index
    is the key for the inner dictionary.
    """

    CDS_dict = {}
    exon_dict = {}
    UTR3_dict = {}
    UTR5_dict = {}
    gene_dict = {}
    ncRNA_dict = {}
    pseudogene_dict = {}

    file_names = ["CDS_merge.bed",
                  "exon_merge.bed",
                  "gene_merge.bed",
                  "ncRNA_merge.bed",
                  "pseudogene_merge.bed",
                  "UTR3_merge.bed",
                  "UTR5_merge.bed"]

    for filename in os.listdir(annotation_dir):

        for filename in file_names:

            filen = os.path.join(annotation_dir, filename)

            with open(filen, 'r') as f:

                for line in f:

                    line_splitted = line.rstrip().split('\t')

                    chromosome = line_splitted[0]
                    start = int(line_splitted[1])
                    end = int(line_splitted[2])

                    if chromosome not in CDS_dict:
                        CDS_dict[chromosome] = {}
                        exon_dict[chromosome] = {}
                        UTR3_dict[chromosome] = {}
                        UTR5_dict[chromosome] = {}
                        gene_dict[chromosome] = {}
                        ncRNA_dict[chromosome] = {}
                        pseudogene_dict[chromosome] = {}

                    if "CDS" in filen:
                        CDS_dict[chromosome][start] = end
                    elif "exon" in filen:
                        exon_dict[chromosome][start] = end
                    elif "3UTR" in filen:
                        UTR3_dict[chromosome][start] = end
                    elif "5UTR" in filen:
                        UTR5_dict[chromosome][start] = end
                    elif "ncRNA" in filen:
                        ncRNA_dict[chromosome][start] = end
                    elif "pseudogene" in filen:
                        pseudogene_dict[chromosome][start] = end
                    elif "gene_merge" in filen:
                        gene_dict[chromosome][start] = end

    return CDS_dict, exon_dict, UTR3_dict, UTR5_dict,\
        gene_dict, ncRNA_dict, pseudogene_dict


def go_through_sample_file(directory,
                           CDS_dict,
                           exon_dict,
                           UTR3_dict,
                           UTR5_dict,
                           gene_dict,
                           ncRNA_dict,
                           pseudogene_dict,
                           tsm_file,
                           background=False,
                           background_tag=None):
    """
    The script go through the given file and
    annotates whether TSM locates in an intergenic
    or in gene regions. By default screens the sample file
    (script searches a filename containing 'tsms' and reads that.)
    If background is True,the script analyses files with
    a name containing 'background'.

    Arguments
    =========

    directory: A directory where files to be analyzed are
               located.
    annotation_dir: A directory where the annotation files are located.
    background: True or False (default: False)


    Results
    =======
    Returns counts for different gene parts.

    """

    CDS_counts = 0
    exon_counts = 0
    intron_counts = 0
    UTR3_counts = 0
    UTR5_counts = 0
    intergenic_counts = 0
    gene_counts = 0
    ncRNA_counts = 0
    pseudogene_counts = 0

    for filename in os.listdir(directory):

        if background is False:

            if tsm_file in filename:

                filen = os.path.join(directory, filename)
            else:
                continue

        elif background is True:

            if background_tag in filename:

                filen = os.path.join(directory, filename)
            else:
                continue

        with open(filen, 'r') as f:

            for line in f:

                line_splitted = line.rstrip().split('\t')
                chromosome = line_splitted[0]
                if 'chr' not in chromosome:
                    chromosome = 'chr' + chromosome
                start = int(line_splitted[1])
                end = int(line_splitted[2])

                found_in_gene = False
                found_in_gene_part = False

                for start_back in gene_dict[chromosome]:

                    end_back = gene_dict[chromosome][start_back]

                    if start_back <= start and end <= end_back:
                        gene_counts += 1
                        found_in_gene = True
                        break
                    elif start <= start_back and start_back < end:

                        gene_counts += 1
                        found_in_gene = True
                        break

                    elif start < end_back and end_back <= end:

                        gene_counts += 1
                        found_in_gene = True
                        break

                    elif start <= start_back and end_back <= end:
                        gene_counts += 1
                        found_in_gene = True
                        break

                if found_in_gene is False:

                    for start_back in ncRNA_dict[chromosome]:

                        end_back = ncRNA_dict[chromosome][start_back]

                        if start_back <= start and end <= end_back:
                            ncRNA_counts += 1
                            found_in_gene = True
                            break

                        elif start <= start_back and start_back < end:

                            ncRNA_counts += 1
                            found_in_gene = True
                            break

                        elif start < end_back and end_back <= end:

                            ncRNA_counts += 1
                            found_in_gene = True
                            break

                        elif start <= start_back and end_back <= end:
                            ncRNA_counts += 1
                            found_in_gene = True
                            break

                if found_in_gene is False:

                    for start_back in pseudogene_dict[chromosome]:

                        end_back = pseudogene_dict[chromosome][start_back]

                        if start_back <= start and end <= end_back:
                            pseudogene_counts += 1
                            found_in_gene = True
                            break

                        elif start <= start_back and start_back < end:

                            pseudogene_counts += 1
                            found_in_gene = True
                            break

                        elif start < end_back and end_back <= end:

                            pseudogene_counts += 1
                            found_in_gene = True
                            break

                        if start <= start_back and end_back <= end:
                            pseudogene_counts += 1
                            found_in_gene = True
                            break

                if found_in_gene is False:

                    intergenic_counts += 1
                    continue

                if found_in_gene is True and found_in_gene_part is False:

                    for start_back in CDS_dict[chromosome]:

                        end_back = CDS_dict[chromosome][start_back]

                        if start_back <= start and end <= end_back:
                            CDS_counts += 1
                            found_in_gene_part = True
                            break
                        elif start <= start_back and start_back < end:

                            CDS_counts += 1
                            found_in_gene_part = True
                            break

                        elif start < end_back and end_back <= end:

                            CDS_counts += 1
                            found_in_gene_part = True
                            break

                        elif start <= start_back and end_back <= end:
                            CDS_counts += 1
                            found_in_gene_part = True
                            break

                if found_in_gene is True and found_in_gene_part is False:

                    for start_back in exon_dict[chromosome]:

                        end_back = exon_dict[chromosome][start_back]

                        if start_back <= start and end <= end_back:
                            exon_counts += 1
                            found_in_gene_part = True
                            break
                        elif start <= start_back and start_back < end:

                            exon_counts += 1
                            found_in_gene_part = True
                            break

                        elif start < end_back and end_back <= end:

                            exon_counts += 1
                            found_in_gene_part = True
                            break

                        if start <= start_back and end_back <= end:
                            exon_counts += 1
                            found_in_gene_part = True
                            break

                if found_in_gene is True and found_in_gene_part is False:

                    for start_back in UTR3_dict[chromosome]:

                        end_back = UTR3_dict[chromosome][start_back]

                        if start_back <= start and end <= end_back:
                            UTR3_counts += 1
                            found_in_gene_part = True

                        elif start <= start_back and start_back < end:

                            UTR3_counts += 1
                            found_in_gene_part = True
                            break

                        elif start < end_back and end_back <= end:

                            UTR3_counts += 1
                            found_in_gene_part = True
                            break

                        if start <= start_back and end_back <= end:
                            UTR3_counts += 1
                            found_in_gene_part = True
                            break

                if found_in_gene is True and found_in_gene_part is False:

                    for start_back in UTR5_dict[chromosome]:

                        end_back = UTR5_dict[chromosome][start_back]

                        if start_back <= start and end <= end_back:
                            UTR5_counts += 1
                            found_in_gene_part = True

                        elif start <= start_back and start_back < end:

                            UTR5_counts += 1
                            found_in_gene_part = True
                            break

                        elif start < end_back and end_back <= end:

                            UTR5_counts += 1
                            found_in_gene_part = True
                            break

                        if start <= start_back and end_back <= end:
                            UTR5_counts += 1
                            found_in_gene_part = True
                            break

                if found_in_gene is True and found_in_gene_part is False:

                    intron_counts += 1
                    continue

    return (CDS_counts, exon_counts, intron_counts, UTR3_counts, UTR5_counts,
            intergenic_counts, gene_counts, ncRNA_counts, pseudogene_counts)


def proportions_z_test(counts_sample_obs, counts_background_obs,
                       sample_trials, background_trials):
    """
    Performs proportions z test.
    """

    from statsmodels.stats.proportion import proportions_ztest

    val = counts_background_obs/background_trials
    prop = proportions_ztest(counts_sample_obs,
                             sample_trials,
                             val,
                             alternative='smaller')

    val = counts_background_obs/background_trials
    prop_larger = proportions_ztest(
            counts_sample_obs,
            sample_trials,
            val,
            alternative='larger')

    return prop, prop_larger


def main():
    """
    Run the enrichment analysis.
    """

    sample_dir = sys.argv[1]
    background_dir = sys.argv[2]
    annotation_dir = sys.argv[3]
    sample_trials = int(sys.argv[4])
    background_trials = int(sys.argv[5])
    tsm_file = sys.argv[6]
    background_file_tag = sys.argv[7]

    (CDS_dict,
     exon_dict,
     UTR3_dict,
     UTR5_dict,
     gene_dict,
     ncRNA_dict,
     pseudogene_dict) = create_annotation_dicts(annotation_dir)

    (CDS_counts,
     exon_counts,
     intron_counts,
     UTR3_counts,
     UTR5_counts,
     intergenic_counts,
     gene_counts,
     ncRNA_counts,
     pseudogene_counts) = go_through_sample_file(sample_dir,
                                                 CDS_dict,
                                                 exon_dict,
                                                 UTR3_dict,
                                                 UTR5_dict,
                                                 gene_dict,
                                                 ncRNA_dict,
                                                 pseudogene_dict,
                                                 tsm_file,
                                                 background=False)

    (CDS_counts_back,
     exon_counts_back,
     intron_counts_back,
     UTR3_counts_back,
     UTR5_counts_back,
     intergenic_counts_back,
     gene_counts_back,
     ncRNA_counts_back,
     pseudogene_counts_back) = go_through_sample_file(
             background_dir,
             CDS_dict,
             exon_dict,
             UTR3_dict,
             UTR5_dict,
             gene_dict,
             ncRNA_dict,
             pseudogene_dict,
             tsm_file,
             background=True,
             background_tag=background_file_tag)

    for item in [(CDS_counts, CDS_counts_back, "CDS"),
                 (exon_counts, exon_counts_back, "exon"),
                 (intron_counts, intron_counts_back, "intron"),
                 (UTR3_counts, UTR3_counts_back, "3UTR"),
                 (UTR5_counts, UTR5_counts_back, "5UTR"),
                 (intergenic_counts, intergenic_counts_back, "intergenic"),
                 (gene_counts, gene_counts_back, "gene"),
                 (ncRNA_counts, ncRNA_counts_back, "ncRNA"),
                 (pseudogene_counts, pseudogene_counts_back, "pseudogene")]:

        prop, prop_larger = proportions_z_test(item[0], item[1],
                                               sample_trials,
                                               background_trials)

        print(item[2], item[0], item[1], Decimal(item[0]/sample_trials),
              Decimal(item[1]/background_trials), 'smaller',
              Decimal(prop[0]), Decimal(prop[1]),
              'larger', Decimal(prop_larger[0]), Decimal(prop_larger[1]))


if __name__ == "__main__":
    main()
