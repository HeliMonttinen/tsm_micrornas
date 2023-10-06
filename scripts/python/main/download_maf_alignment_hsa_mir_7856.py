""""
A modified script from the download_maf_alignment function
to create a fig from hsa-mir-7856 as a few primates contain
a long insertion on left-hand-side, which we want to include to the
fig.

"""
from Bio import AlignIO
from collections import defaultdict
import sys


def locus_ids(data_file):
    """
    Collects unique identifiers for each
    given locus, because multiple identifiers maybe
    associated the same locus.

    Argument
    ========
    :data_file: A data file containing the
                informationtin on identifier and
                locus separated by a tab. An identifier
                is located in the first column and the
                locus in the last.
    Output
    ======
    Identifiers: A set containing unique identifiers.
    """

    locus_dict = {}
    identifiers = {}

    with open(data_file, 'r') as f:
        for line in f:
            if not 'hsa-mir-7856' in line:
                continue
            line_splitted = line.rstrip().split()
            EGN_identifier = line_splitted[0]
            locus = line_splitted[-1].rstrip()
            chromosome = locus.split(':')[0]
            start = int(locus.split(':')[1].split('-')[0])-290
            end = int(locus.split(':')[1].split('-')[1])
            reg = set(range(start, end))
            found_same_seq = False

            for loc in locus_dict:
                chr_set = loc.split(':')[0]
                if chr_set == chromosome:
                    start_set = int(loc.split(':')[1].split('-')[0])-290
                    end_set = int(loc.split(':')[1].split('-')[1])
                    set_reg = set(range(start_set, end_set))
                    intersec = reg.intersection(set_reg)
                    if len(intersec) > 0:
                        print("bad", EGN_identifier)
                        found_same_seq = True
                        break

            if found_same_seq is False:
                locus_dict[locus] = EGN_identifier

    for locs in locus_dict:
        chr_n = str(locs.split(':')[0])
        start = int((locs.split(':')[-1]).split('-')[0])-290
        end = int(locs.split('-')[-1])
        identifiers[locus_dict[locs]] = [chr_n, start, end]

    return identifiers


def spes_names(species_file):
    """
    Returns species name dict

    Argument: species_file

    Returns:
    A dictionary, in which genome identifier is the key
    and species name a value

    """
    species_dict = {}
    with open(species_file, 'r') as f:

        for line in f:
            if '=' in line:
                continue
            if len(line.rstrip()) == 0:
                continue
            species_list = line.rstrip().split()[1:]
            species = ""
            upper_letter = False
            for spes in species_list:
                if spes[0].isupper() and upper_letter is False:
                    upper_letter = True
                    species = species + ' ' + spes
                elif upper_letter is True and not spes[0].isupper():
                    species = species + ' ' + spes
                elif spes[0].isupper() and upper_letter is True:

                    identifier = line.split('/')[-1].split()[0]
                    species_dict[identifier] = species
                    break

    return species_dict


def extract_areas(start, end, identifier, chr_file,
                  tar_dir, maf_dir, species_dict, seq_info_dir):
    """
    Extracts the coordinates and sequences for different species
    from .maf alignments.


    Arguments
    =========

    start: A staring index of a microRNA locus
    end: The last index of a microRNA locus
    identifier: identifier of a microRNA
    chr_file: .maf alignment file
    tar_dir: A path to the target dictionary, where sequences
             and species names are saved.
    maf_dir: A path to a directory, where the extracted
             .maf alignment is saved.
    species_dict: A dictionary, where genome assembly identifiers
                  are keys and species names as values.
                  The output of a function 'spes_names'.
    seq_info_dir: A path to a directory, where information
                  on the sequences is saved. The file is tab-separated
                  and columns contain following information:
                  identifies, species_name, fragment_size,
                  a strand of the starting index,
                  the strand of the last index, a starting index,
                  the last index

    Output
    ======
    files in tar_dir, maf_dir and seq_info_dir

    """
    index_file = chr_file.split('.')[0] + '.mafindex'
    target = "hg38." + (chr_file.split('.')[0]).split('/')[-1]
    idx = AlignIO.MafIO.MafIndex(index_file, chr_file, target)
    results = idx.search([start], [end])

    block = 0
    info = {}
    block_count = defaultdict(dict)
    len_dict = {}
    found_spes = []
    seq_info_dict = defaultdict(dict)

    identifier = str(identifier)
    for result in results:

        with open(maf_dir + identifier + '.fas', 'a+') as f1:
            f1.write(str(result))

        block += 1
        block_count[block] = dict()
        for line in result:
            if 'Alignment' in line:
                block_count += 1
            if 'Alignment' not in line:
                header = str(line).rstrip().split()[1].split('.')[0]
                if header not in species_dict:
                    continue
                if header not in found_spes:
                    found_spes.append(header)
                    info_dict = {"fragment_size": 0,
                                 "start_strand": 'None',
                                 "end_strand": 'None'}

                    seq_info_dict[header] = info_dict

                block_count[block][header] = str(line.seq)
                len_dict[block] = len(line.seq)
                chr_name = '.'.join(str(line.name).rstrip().split('.')[1:])
                fragment_size = str(line).split("srcSize=")[-1].split()[0]
                seq_info_dict[header]["fragment_size"] = fragment_size
                strand = str(line).split("strand=")[-1].split()[0]

                if '-' in strand:
                    seq_info_dict[header]["end_strand"] = "-"
                    if seq_info_dict[header]["start_strand"] ==\
                            "None":
                        seq_info_dict[header]["start_strand"] =\
                                "-"
                else:
                    seq_info_dict[header]["end_strand"] = "+"
                    if seq_info_dict[header]["start_strand"] == "None":
                        seq_info_dict[header]["start_strand"] = "+"

                start = str(line).split("start=")[-1].split()[0]

                if header not in info:
                    end = int(start) +\
                            int((str(line).split("size=")[-1]).split()[0])
                    info_line = chr_name + ':' + start + '_' + str(end)
                    info[header] = info_line
                else:
                    end = int(start) +\
                            int((str(line).split("size=")[-1]).split()[0])
                    info_line = info[header]
                    info_line_prel = '_'.join(info_line.split('_')[:-1])
                    info_line = info_line_prel + '_' + str(end)
                    info[header] = info_line

    final = {}

    for ids in block_count:
        for head in found_spes:

            if head in block_count[ids]:
                if head in final:

                    final[head] = final[head] + block_count[ids][head]
                else:
                    final[head] = block_count[ids][head]

            else:
                if head in final:

                    final[head] = final[head] + '-'*len_dict[ids]
                else:
                    final[head] = '-'*len_dict[ids]

    final2 = {}

    for item in final:
        max_ids = 0

        for item2 in final:

            if item == item2:
                continue

            ids = 0
            ids2 = 0
            align_len = 0

            for a,b in zip(final[item].upper(),
                           final[item2].upper()):

                if a == b and a != '-':
                    ids += 1
                if a != '-' or b != '-':
                    align_len += 1

            ids2 = ids/align_len
            if float(ids2) > float(max_ids):
                max_ids = ids2

        if max_ids > 0.40 and len(final[item].replace(
                '-', '')) > 40:
            final2[item] = final[item]
        elif 'hg38' in item:
            final2[item] = final[item]

    with open(seq_info_dir + identifier + '.txt', '+a') as f1:

        for spes in seq_info_dict:

            if str(seq_info_dict[spes]["start_strand"]) == '+':
                start_index = info[spes].replace(
                        ':', '_').split('_')[-2]
            else:
                frag = int(seq_info_dict[spes]["fragment_size"])
                start_index = frag - int(info[spes].replace(
                    ':', '_').split('_')[-2])

            if str(seq_info_dict[spes]["end_strand"]) == '+':
                end_index = info[spes].replace(
                        ':', '_').split('_')[-1]

            else:
                frag = int(seq_info_dict[spes]["fragment_size"])
                end_index = frag - int(info[spes].replace(
                    ':', '_').split('_')[-1])

            line_text = identifier + '\t' + species_dict[spes] + '\t' +\
                str(seq_info_dict[spes]["fragment_size"]) + '\t' +\
                str(seq_info_dict[spes]["start_strand"]) + '\t' +\
                str(seq_info_dict[spes]["end_strand"]) + '\t' +\
                str(start_index) + '\t' +\
                str(end_index) + '\n'

            f1.write(line_text)

            with open(tar_dir + identifier + '.fas', 'a+') as f2:

                if spes in final2:
                    seq = (final2[spes]).upper() + '\n'
                    new_header = '>' + (species_dict[spes]).replace(
                            ' ', '_') + '_' + info[spes].replace(
                                    ':', '_') + '\n'
                    f2.write(new_header)
                    f2.write(seq)
    return


def main():
    """
    Retrieves corresponding coordinates from .maf alignments
    for all microRNA coordinates in the given data file and saves them
    into separate files.

    Arguments
    ==========
    :data_file: A tab-separated microRNAfile containing
                identifiers and coordinates. Identifiers are
                in the first column, coordinates in the last.
    :species_file: A .maf alignment associated file containing
                   species and their assemblies.
    :chr_dir: A path to a directory, where .maf alignment files are
              located.
    :target_dir: A path to a directory, where extracted sequences
                 sequences are saved.
    :maf_dir: A path to a directory, where extract .maf alignments
              are saved.
    seq_info_dir: A path to directory, where information on
                  the extracted .maf alignment sequences are saved
                  (identifier, species, seq length, strands of
                   starting and ending coordinates and coordinates).

    """

    data_file = sys.argv[1]
    species_file = sys.argv[2]
    chr_dir = sys.argv[3]
    target_dir = sys.argv[4]
    maf_dir = sys.argv[5]
    seq_info_dir = sys.argv[6]

    identifiers = locus_ids(data_file)
    species_names = spes_names(species_file)

    for identifier in identifiers:
        if len(identifiers[identifier][0]) > 5:
            continue

        chr_file = chr_dir + 'chr' + identifiers[identifier][0] + '.maf'

        start = int(identifiers[identifier][1])
        end = int(identifiers[identifier][2])
        extract_areas(
                start,
                end,
                identifier,
                chr_file,
                target_dir,
                maf_dir,
                species_names,
                seq_info_dir)


if __name__ == "__main__":
    main()
