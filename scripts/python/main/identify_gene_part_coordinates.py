from collections import (defaultdict,
                         OrderedDict)
import sys

def get_coordinates(line):
    """
    Parses coordinates from a line.
    """

    line_splitted = line.split()
    chrom = line_splitted[0]
    start = int(line_splitted[3])
    stop = int(line_splitted[4])
    if "Parent=gene:" in line:
        parent = line_splitted[8].split('Parent=gene:')[-1].split(';')[0]
    else:
        parent = None
    return parent, chrom, start, stop


def read_gff_file(gff_file):
    """
    """

    chr_found = False
    pargene = False

    gene_info_dict_exon = defaultdict(dict)
    gene_info_dict_CDS = defaultdict(dict)
    gene_info_dict_UTR3 = defaultdict(dict)
    gene_info_dict_UTR5 = defaultdict(dict)
    parent_type_dict_ncRNA = defaultdict(dict)
    parent_type_dict_gene = defaultdict(dict)
    parent_type_dict_pseudogene = defaultdict(dict)

    with open(gff_file, 'r') as f:

        gene_found = False
        for line in f:

            if 'ID=gene:' in line:
                pargene = True
                parent_name = line.split()[-7].split(';')[0].lstrip('ID=gene:')
                if 'chr' + line.split()[0]  not  in parent_type_dict_ncRNA:
                    parent_type_dict_ncRNA['chr' + line.split()[0]] = OrderedDict()
                    parent_type_dict_gene['chr' + line.split()[0]] = OrderedDict()
                    parent_type_dict_pseudogene['chr' + line.split()[0]] = OrderedDict()

                if line.split()[2] == "ncRNA_gene":

                    parent_type_dict_ncRNA[
                            'chr' + line.split()[0]][int(line.split()[3]),
                            int(line.split()[4])] = line.split()[2]

                elif line.split()[2] == "gene":

                    parent_type_dict_gene[
                            'chr' + line.split()[0]][int(line.split()[3]),
                            int(line.split()[4])] = line.split()[2]

                elif line.split()[2] == "pseudogene":

                    parent_type_dict_pseudogene[
                            'chr' + line.split()[0]][int(line.split()[3]),
                            int(line.split()[4])] = line.split()[2]

            if ('###' not in line) and len(line.split()) < 3:
                continue

            if len(line.split()) > 7 and ('exon' in line.split()[2] or\
                    'UTR' in line.split()[2] or 'CDS' in line.split()):

                gen_parent, gen_chrom, gen_start, gen_stop = get_coordinates(line)
                element = line.split()[2]

                gen_chrom = 'chr' + gen_chrom
                result_list = [element, gen_parent]
                if gen_chrom not in gene_info_dict_exon:
                    gene_info_dict_exon[gen_chrom] = OrderedDict()
                    gene_info_dict_CDS[gen_chrom] = OrderedDict()
                    gene_info_dict_UTR3[gen_chrom] = OrderedDict()
                    gene_info_dict_UTR5[gen_chrom] = OrderedDict()

                if element == 'exon':

                    gene_info_dict_exon[gen_chrom][gen_start, gen_stop] = result_list

                if element == "CDS":

                    gene_info_dict_CDS[gen_chrom][gen_start, gen_stop] = result_list

                if element == "three_prime_UTR":

                    gene_info_dict_UTR3[gen_chrom][gen_start, gen_stop] = result_list

                if element == "five_prime_UTR":

                    gene_info_dict_UTR5[gen_chrom][gen_start, gen_stop] = result_list

    return parent_type_dict_ncRNA, parent_type_dict_gene,\
            parent_type_dict_pseudogene,\
            gene_info_dict_exon, gene_info_dict_CDS,\
            gene_info_dict_UTR3, gene_info_dict_UTR5


def write_bed_files(parent_type_dict_ncRNA,
                    parent_type_dict_gene,
                    parent_type_dict_pseudogene,
                    gene_info_dict_exon,
                    gene_info_dict_CDS,
                    gene_info_dict_UTR3,
                    gene_info_dict_UTR5):
    """
    Write separate bed files for each genetic components.

    Argument
    ========
    gene_info_dict: A nested dictionary containing information
                    gene elements
    parent_type_dict: A nested dictionary containing information
                      on gene type (gene or ncRNA gene)
    """

    order = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
             'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2',
             'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
             'chr8', 'chr9', 'chrX', 'chrY']

    for chrom in order:

        help_dict = defaultdict(list)
        for coord in parent_type_dict_ncRNA[chrom]:


            help_dict[coord[0]].append(coord[1])

        keys = sorted(help_dict.keys())

        for key in keys:

            with open('results/genome_annotations/' + 'ncRNA' + '.bed', 'a+') as f:

                for a in help_dict[key]:

                    f.write(chrom + '\t' + str(key) + '\t' + str(a) + '\n')

        help_dict = defaultdict(list)
        for coord in parent_type_dict_gene[chrom]:


            help_dict[coord[0]].append(coord[1])

        keys = sorted(help_dict.keys())

        for key in keys:

            with open('results/genome_annotations/' + 'gene' + '.bed', 'a+') as f:
                for a in help_dict[key]:

                    f.write(chrom + '\t' + str(key) + '\t' + str(a) + '\n')

        help_dict = defaultdict(list)
        for coord in parent_type_dict_pseudogene[chrom]:

            help_dict[coord[0]].append(coord[1])

        keys = sorted(help_dict.keys())

        for key in keys:

            with open('results/genome_annotations/' + 'pseudogene' + '.bed', 'a+') as f:

                for a in help_dict[key]:

                    f.write(chrom + '\t' + str(key) + '\t' + str(a) + '\n')

        help_dict = defaultdict(list)
        for coord in gene_info_dict_exon[chrom]:

            help_dict[coord[0]].append(coord[1])

        keys = sorted(help_dict.keys())

        for key in keys:

            with open('results/genome_annotations/' + 'exon' + '.bed', 'a+') as f:

                for a in help_dict[key]:

                    f.write(chrom + '\t' + str(key) + '\t' + str(a) + '\n')

        help_dict = defaultdict(list)
        for coord in gene_info_dict_CDS[chrom]:
            
            help_dict[coord[0]].append(coord[1])

        keys = sorted(help_dict.keys())

        for key in keys:

            with open('results/genome_annotations/' + 'CDS' + '.bed', 'a+') as f:

                for a in help_dict[key]:

                    f.write(chrom + '\t' + str(key) + '\t' + str(a) + '\n')

        help_dict = defaultdict(list)
        for coord in gene_info_dict_UTR3[chrom]:

            help_dict[coord[0]].append(coord[1])

        keys = sorted(help_dict.keys())

        for key in keys:

            with open('results/genome_annotations/' + 'UTR3' + '.bed', 'a+') as f:

                for a in help_dict[key]:

                    f.write(chrom + '\t' + str(key) + '\t' + str(a) + '\n')

        help_dict = defaultdict(list)
        for coord in gene_info_dict_UTR5[chrom]:

            help_dict[coord[0]].append(coord[1])

        keys = sorted(help_dict.keys())

        for key in keys:

            with open('results/genome_annotations/' + 'UTR5' + '.bed', 'a+') as f:

                for a in help_dict[key]:

                    f.write(chrom + '\t' + str(key) + '\t' + str(a) + '\n')

def main():
    """
    Run the script
    """


    gff3_file = sys.argv[1]


    (parent_type_dict_ncRNA,
     parent_type_dict_gene,
     parent_type_dict_pseudogene,
     gene_info_dict_exon,
     gene_info_dict_CDS,
     gene_info_dict_UTR3,
     gene_info_dict_UTR5) = read_gff_file(gff3_file)

    write_bed_files(parent_type_dict_ncRNA,
                    parent_type_dict_gene,
                    parent_type_dict_pseudogene,
                    gene_info_dict_exon,
                    gene_info_dict_CDS,
                    gene_info_dict_UTR3,
                    gene_info_dict_UTR5)


if __name__ == "__main__":
    main()
