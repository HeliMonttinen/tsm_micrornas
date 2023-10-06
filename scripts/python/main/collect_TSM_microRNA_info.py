"""
This script collects information on microRNAs,
which locus contain a TSM longer than 15 bases causing
at least five mismatches. Identifies the location of the
mir in relation to other genes (intron, exon, protein-coding gene,
non-coding gene). Finally writes all the information to a matrix.
"""

from collections import defaultdict
import os
import sys
from pymongo import MongoClient


dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def parse_dat_file(dat_file):
    """
    Collects microRNA information (miRbase accession and
    the mature microRNa locus coordinates from a .dat file.

    Arguments
    =========
    dat_file: A .dat file downloaded from miRBase

    Returns
    ========
    mir_products: A dictionary where an identifier is a key
                  and value is a dictionary with information
                  on start and stop indexes.
    """

    mir_products = defaultdict(dict)

    same = False

    with open(dat_file, 'r') as f:
        for line in f:

            if '//' in line:
                if same is True:
                    mir_products[identifier] = MIRNA_acc
                MIRNA_acc = dict()
                same = False
                MIRNA_acc = dict()

            if 'ID' in line and 'hsa' in line:
                same = True

                identifier = line.rstrip().split()[1]

            if 'FT' in line and 'miRNA' in line and same is True:

                start = line.split()[-1].split('..')[0]
                stop = line.rstrip().split()[-1].split('..')[1]

            if 'FT' in line and 'accession' in line and same is True:
                accession = line.rstrip().split('=')[-1]
                MIRNA_acc[accession] = {"start": start,
                                        "stop": stop}

    return mir_products


def connect_EGNS_to_mir(EGNS_file, mir_products, database, output,
                        genomic_context):
    """
    Connects information on microRNA Ensemble identifier,
    microRNA name and genomic context. Writes the information
    into a .csv file

    Arguments
    =========

    EGNS_file: A path to file containing information on
               Ensemble identifier, microRNA name and locus
    mir_products: A dictionary containing information on mature microRNAs
                  and their coordinates.
    database: The name of the mongo database where TSM-associated
               microRNAs are saved.
    output: Where the output is saved
    Genomic_context: Information parsed from a .gff file and
                     retrieved GO terms for the host gene
    """

    locus_dict = defaultdict(list)
    id_dict = defaultdict(list)

    with open(EGNS_file, 'r') as f:
        for line in f:

            line_splitted = line.rstrip().split()
            if len(line_splitted) == 4:

                locus_dict[line_splitted[-1]].append(line_splitted[0])
                id_dict[line_splitted[0]].append(line_splitted[2])

    client = MongoClient()
    db = client[database]
    collection = db['TSMs']

    query = {"mismatches": {"$gt": 5},
             "TSM_length": {"$gt": 15},
             "quality": True}

    for document in collection.find(query):

        hit_dict = {}

        for key, value in document.items():

            hit_dict[key] = value

        location = hit_dict["human_location"]

        accession = ""
        products = ""

        identifiers_same_loc = ', '.join(locus_dict[location])
        used = set()
        accession_set = set()

        for other_loc in locus_dict[location]:

            for mir_loc in id_dict[other_loc]:
                accession_set.add(mir_loc)

                for mir_pro in mir_products[mir_loc]:

                    if mir_pro in used:
                        continue
                    used.add(mir_pro)

                    if len(products) == 0:
                        products = mir_pro.lstrip('"').rstrip('"') + ':' +\
                            '(' + mir_products[mir_loc][mir_pro]["start"] +\
                            '-' + mir_products[mir_loc][mir_pro]["stop"] + ')'
                    else:
                        products = products + ', ' +\
                            mir_pro.lstrip('"').rstrip('"') + ':' +\
                            '(' + mir_products[mir_loc][mir_pro]["start"] +\
                            '-' + mir_products[mir_loc][mir_pro]["stop"] + ') '

        accession = ', '.join(list(accession_set))

        species = ', '.join(hit_dict["taxonomic_group"])

        host_info = genomic_context[hit_dict["identifier"]]

        if len(host_info) == 0:

            host_info = {'mRNA': "None",
                         'lncRNA': "None",
                         'microRNA': "None",
                         'parent_gene_identifier_mRNA': set(),
                         'parent_gene_identifier_lncRNA': set(),
                         'exon': "None",
                         'BP_ids': "None",
                         'BP_terms': "None",
                         'CC_ids': "None",
                         'CC_terms': "None",
                         'MF_ids': "None",
                         'MF_terms': "None",
                         'host_identifier': "None",
                         'host_description': "None",
                         'UTR': "None"}

        output_text = hit_dict["identifier"] + '\t' +\
            hit_dict["human_location"] + '\t' +\
            str(hit_dict["TSM_length"]) + '\t' +\
            identifiers_same_loc + '\t' +\
            accession + '\t' + products + '\t' +\
            host_info["microRNA"] + '\t' + host_info["lncRNA"] + '\t' +\
            host_info["mRNA"] + '\t' + host_info["UTR"] +\
            '\t' + host_info["exon"] + '\t' +\
            host_info["host_identifier"] + '\t' +\
            host_info["host_description"] + '\t' +\
            host_info["BP_ids"] + '\t' + host_info["BP_terms"] + '\t' +\
            host_info["CC_ids"] + '\t' + host_info["CC_terms"] + '\t' +\
            host_info["MF_ids"] + '\t' + host_info["MF_terms"] +\
            '\t' + species + '\n'

        with open(output, 'a+') as f1:
            f1.write(output_text)


def genomic_context(gff_info):
    """
    Collects information on the genomic context of a microRNA.
    If the microRNA is located within another gene,
    it retrieves GO terms for the host gene.

    Argeuments
    ===========
    gff_info: A path to a parsed gff file where microRNA id and
              locus are in twp first columns.

    Returns
    =======

    genomic_context_dict: A dictionary where microRNA
                          identifier is a key and combined information
                          a value in a string format.

    """

    from go_terms import get_go_terms
    genomic_context_dict = defaultdict(dict)

    with open(gff_info, 'r') as f1:

        for line in f1:

            line_splitted = line.rstrip().split('\t')
            for item in line_splitted:
                if len(item) == 0:
                    line_splitted.remove("")

            if line_splitted[0] not in genomic_context_dict:

                info = {'gene': "None",
                        'ncRNA_gene': "None",
                        'UTR': "None",
                        'mRNA': "None",
                        'lncRNA': "None",
                        'microRNA': "",
                        'exon': "None",
                        'BP_ids': "",
                        'BP_terms': "",
                        'CC_ids': "",
                        'CC_terms': "",
                        'MF_ids': "",
                        'MF_terms': "",
                        'host_identifier': "None",
                        'host_description': "None",
                        'UTR': "None"}

            if 'havana' in line_splitted[4] or\
                    'ensembl' in line_splitted[4]:

                if 'lnc_RNA' in line_splitted[5]:

                    if info["lncRNA"] == "None":
                        info["lncRNA"] = line_splitted[9]
                    elif line_splitted[9] not in info["lncRNA"]:
                        info["lncRNA"] = info["lncRNA"] + '/' +\
                                line_splitted[9]

                elif 'mRNA' in line_splitted[5]:
                    if info["mRNA"] == "None":
                        info["mRNA"] = line_splitted[9].rstrip().lstrip()
                    elif line_splitted[9] not in info["mRNA"]:
                        info["mRNA"] = info["mRNA"] + '/' + line_splitted[9]

                elif "exon" in line_splitted[5]:
                    if info["exon"] == "None":
                        info["exon"] = line_splitted[9]
                    elif info["exon"] != line_splitted[9]:
                        info["exon"] = info["exon"] + '/' + line_splitted[9]

                elif "gene" == line_splitted[5] or\
                        "ncRNA_gene" == line_splitted[5]:
                    host_identifier_sp =\
                            line_splitted[11].split(';')[0].split("ID=gene:")[-1]
                    if 'description' in line_splitted[11].split(';')[3]:
                        host_description_sp =\
                            line_splitted[11].split(
                                ';')[3].split('[')[0].split('description=')[-1]
                    else:
                        host_description_sp =\
                                line_splitted[11].split(
                                    ';')[2].split('=')[-1].split('[')[0]

                    if info["host_identifier"] != "None":
                        info["host_identifier"] =\
                            info["host_identifier"] + '/' + host_identifier_sp
                        info["host_description"] =\
                            info["host_description"] + '/' +\
                            host_description_sp
                    else:
                        info["host_identifier"] = host_identifier_sp
                        info["host_description"] = host_description_sp

                    go_ids, go_terms, _ = get_go_terms(host_identifier_sp)

                    if len(go_ids["BP"]) > 0:
                        info["BP_ids"] = info["BP_ids"] +\
                            ', '.join(list(set(go_ids["BP"])))
                        info["BP_terms"] = info["BP_terms"] +\
                            ', '.join(list(set(go_terms["BP"])))
                    if len(go_ids["CC"]) > 0:
                        info["CC_ids"] = info["CC_ids"] +\
                            ', '.join(list(set(go_ids["CC"])))
                        info["CC_terms"] = info["CC_terms"] +\
                            ', '.join(list(set(go_terms["CC"])))
                    if len(go_ids["MF"]):
                        info["MF_ids"] = info["MF_ids"] +\
                            ', '.join(list(set(go_ids["MF"])))
                        info["MF_terms"] = info["MF_terms"] +\
                            ', '.join(list(set(go_terms["MF"])))

                    if len(info["BP_ids"]) == 0:
                        info["BP_ids"] = "None"
                        info["BP_terms"] = "None"
                    if len(info["CC_ids"]) == 0:
                        info["CC_ids"] = "None"
                        info["CC_terms"] = "None"
                    if len(info["MF_ids"]) == 0:
                        info["MF_ids"] = "None"
                        info["MF_terms"] = "None"

                elif "UTR" in line_splitted[5]:
                    info["UTR"] = line_splitted[9]

            elif 'mirbase' in line_splitted[4]:

                if len(info["microRNA"]) == 0:
                    info["microRNA"] = line_splitted[9]
                elif line_splitted[9] not in info["microRNA"]:
                    info["microRNA"] = info["microRNA"] +\
                            '/' + line_splitted[9]

            genomic_context_dict[line_splitted[0]] = info

    return genomic_context_dict


def main():
    """
    Runs the code.

    Arguments
    =========

    dat_file: a RNA .dat file downloaded from miRbase
    EGNS_file: A -csv file in which Ensemble id, Mirbase access,
               mirname, and coordinates
    output filename: Where the file is saved
    database: A database, in which TSM microRNAs are saved
    gff_info: .A file path to -.gff file
    """

    dat_file = sys.argv[1]
    EGNS_file = sys.argv[2]
    output = sys.argv[3]
    database = sys.argv[4]
    gff_info = sys.argv[5]

    mir_products = parse_dat_file(dat_file)

    genomic_context_dict = genomic_context(gff_info)

    connect_EGNS_to_mir(EGNS_file,
                        mir_products,
                        database,
                        output,
                        genomic_context_dict)


if __name__ == "__main__":
    main()
