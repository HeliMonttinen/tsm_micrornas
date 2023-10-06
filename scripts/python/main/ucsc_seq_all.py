"""
Extracts the full sequence from genomes based on the
given genomic coordinates
"""
from collections import defaultdict
import os
import subprocess
import sys
from Bio.Seq import reverse_complement

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def parse_alignment_file(data):
    """
    Yield the location of the sequence in a genome.

    Arguments
    ==========
    data: The file  containing
          species names and assembly identifiers.

    Yield
    =====
    accession:Assembly identifier
    species_name: A scientific species name
    The sequence location
    """
    with open(data, 'r') as f:

        for line in f:
            if '>' in line:

                accession = data.split('/')[-1].split('.')[0]
                line_splitted = line.lstrip('>').rstrip().split('_')
                if 'chr' in line_splitted[-5] and len(line_splitted[-5]) < 7:
                    species_name = '_'.join(line_splitted[:-5])
                    chromosome = '_'.join(line_splitted[-5:-2])

                elif line_splitted[-4] == 'scaffold' or\
                        ('chr' in line_splitted[-4] and
                         len(line_splitted[-4]) < 7):
                    species_name = '_'.join(line_splitted[:-4])
                    chromosome = '_'.join(line_splitted[-4:-2])
                else:
                    species_name = '_'.join(line_splitted[:-3])
                    chromosome = line_splitted[-3]
                start = line_splitted[-2]
                end = line_splitted[-1]

                yield (accession, species_name,
                       chromosome + ':' + str(start) + '-' + str(end))


def parse_start_end_sites(locus_indexes_file):
    """
    Identify start and stop indexes for each sequence.
    Ignore if the first and last .maf synteny blocks are
    on different strands.

    Arguments
    =========

    locus_indexes_file: The .synteny block file.
    """

    start_stop_dict = defaultdict(list)
    reverse_dict = {}

    with open(locus_indexes_file, 'r') as f:

        for line in f:
            line_splitted = line.split('\t')
            if line_splitted[-3] != line_splitted[-4]:
                continue
            species_name = line_splitted[1].lstrip('_').lstrip()
            number1 = (line_splitted[-2])
            number2 = (line_splitted[-1].rstrip())
            start_stop_dict[species_name].append(number1)
            start_stop_dict[species_name].append(number2)
            start_stop_dict[species_name] =\
                list(set(start_stop_dict[species_name]))
            start_stop_dict[species_name].sort()
            reverse_dict[species_name] = line_splitted[-3]

    return start_stop_dict, reverse_dict


def get_assembly_id(assembly_ids):
    """
    Get assemble identifiers for species and
    create a dictionary.

    Argument
    =========

    assembly_ids: The species_mammals.txt file
    """

    assembly_id_dict = {}

    with open(assembly_ids, 'r') as f:

        for line in f:

            if len(line) > 2 and '=' not in line:

                species_list = line.rstrip().split()[1:]
                species = ""
                upper_letter = False
                for spes in species_list:
                    if spes[0].isupper() and upper_letter is False:
                        upper_letter = True
                        species = species + '_' + spes
                    elif upper_letter is True and not spes[0].isupper():
                        species = species + '_' + spes
                    elif spes[0].isupper() and upper_letter is True:
                        identifier = line.split('/')[-1].split()[0]
                        species = species.lstrip('_')
                        assembly_id_dict[species] = identifier
                        break

    return assembly_id_dict


def main():
    """
        Screens the data directory for loci sequences in a fasta format.
    Retrieves the full length sequences
    from the the genome assemblies based on the species name,
    chromosome and starting and ending indexes.
    The genome assembly identifiers are parsed from
    the annotation file downloaded from the UCSC database.
    Sequences, whose starting and ending indexes are on the opposite strands
    are ignored.

    Arguments
    =========

    datadir: The directory containing sequencefiles and
             the files are named as Ensembl_identifier_chunknumber.xxx
    assembly_file: Edited
      http://hgdownload.cse.ucsc.edu/goldenPath/hg38/multiz100way/README.txt
      file. Contains information on the species names and the corresponding
      assembly identifier used in the UCSC .maf alignments.

    info_dir: A full path to a directory containing information files on
              the sequences extracted from the .maf alignments.
              Files are tab-separated and contain information on:
              species, chromosome, starting and ending indexes,
              and strand info.


    Output
    ======

    Full length sequences saved in a fasta format into the
    results/outputdir/ directory. Each file corresponds
    one gene chunk and they are named with Ensembl_identifier_chunknumber.
    Sequence headers contain a scientific species name,
    chromosome, and starting and stopping indexes.

    """

    assembly_file = sys.argv[1]
    info_dir = sys.argv[2]
    filename = sys.argv[3]
    outputdir = sys.argv[4]

    assembly_id_dict = get_assembly_id(assembly_file)

    filepath = filename
    identifier = filename.split('/')[-1].rstrip('.txt').rstrip('.fas')

    start_stop_dict, reverse_dict = parse_start_end_sites(
            info_dir + '/' + identifier + '.txt')

    for accession, species_name, location in parse_alignment_file(filepath):
        species_name = species_name.lstrip('_')
        species_name2 = species_name.replace('_', ' ')

        try:
            start = str(start_stop_dict[species_name2][0])
            end = str(start_stop_dict[species_name2][1])
        except:
            continue

        chrom = location.split(':')[0]
        if species_name not in assembly_id_dict:
            continue
        assembly = assembly_id_dict[species_name]

        if assembly in ['capHir1', 'monDom5', 'ficAlb2']:
            continue
        if 'random' in chrom:

            chrom = chrom.split('_')[1]

            with open('data/genomes/' +
                      assembly + '.fasta.gz.fai', 'r') as f3:
                for line in f3:

                    line_split = line.rstrip().split('\t')[0]
                    if chrom in line_split:
                        chrom = line_split

        elif 'chr' not in chrom and 'scaffold' not in chrom:
            with open('data/genomes/' +
                      assembly + '.fasta.gz.fai', 'r') as f3:
                for line in f3:

                    line_split = line.rstrip().split('\t')[0]
                    if chrom in line_split:
                        chrom = line_split

        try:

            p = subprocess.check_output([
                "samtools",
                "faidx",
                "data/genomes/{assembly}.fasta.gz".format(
                    assembly=assembly),
                "{chrom}:{start}-{end}".format(
                    chrom=chrom,
                    start=start,
                    end=end)])
        except:
            print(accession, species_name, assembly, chrom, start, end)
            continue

        seq_list = str(p, 'utf-8').rstrip().split('\n')

        seq = ''.join(seq_list[1:]).upper()
        if len(seq) < 200:
            continue

        with open(outputdir + "/" + identifier, 'a+') as f:

            f.write(">" + species_name + '_' + chrom + '_' + start +
                    '_' + end + '\n')

            if reverse_dict[species_name2] == '-':

                seq = reverse_complement(seq)

            f.write(seq + '\n')


if __name__ == "__main__":
    main()
