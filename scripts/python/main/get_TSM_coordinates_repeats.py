"""
This script extracts the TSM coordinates
from a database and save them into a bed formatted file
if TSMs target or source sites do not overlap with
a repeat more than 50%.
"""

from collections import (defaultdict,
                         OrderedDict)
import os
import sys

from Bio import SeqIO
from pymongo import MongoClient
import subprocess

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def run_dustmasker(seq_file,
                   alignment_seq,
                   source_start,
                   source_end,
                   target_start,
                   target_end):
    """
    Run dustmasker for the given sequence and identify
    if TSM regions overlap with a repeat more than 50 of the
    length of a TSM source or target site. If not, return
    False, otherwise not.

    Arguments
    ==========
    seq_file: A file containing non-aligned full length sequence.
              The script searches only Homo sapiens sequence and runs
              Dustmasker for that.
    alignment_seq: Aligned sequence of the corresponding Homo sapiens sequence.
    source_start: A starting index of a source site (aligned sequence)
    source_end: An last index of a source site (aligned sequence)
    target_start: A starting index of a target site (aligned sequence)
    target_end: A last index of a target site (aligned sequence)

    Returns
    ========

    repeat_found: True, if a TSM target or source site overlaps more then 50%
                   with a repeat.
    """
    from RNA_alignments import indexes_in_alignment

    p = subprocess.check_output(["dustmasker",
                                 "-in",
                                 "{seq_file}".format(seq_file=seq_file)])

    line_list = str(p, 'utf-8').split('\n')
    identifier_set = set()
    Homo_sapiens = False
    for line in line_list:
        if '>' in line and 'Homo' in line:
            Homo_sapiens = True
        elif '>' not in line and Homo_sapiens is True:
            start_stop = line.rstrip().split('-')
            start = int(start_stop[0])
            stop = int(start_stop[1])
            start_stop = indexes_in_alignment(
                    [start,
                     stop], alignment_seq)
            identifier_set.add((start_stop[0], start_stop[1]))
        elif '>' in line and Homo_sapiens is True:
            break

    repeat_found = False
    source_site = len(alignment_seq[
        source_start:source_end].replace('-', ''))
    target_site = len(alignment_seq[
        target_start:target_end].replace('-', ''))
    if source_site == 0:
        return True
    elif target_site == 0:
        return True

    for reps in identifier_set:

        if reps[0] <= source_start <= reps[1]:
            len1 = len(alignment_seq[
                source_start:reps[1]].replace('-', ''))
            if len1/source_site > 0.5:

                repeat_found = True
                break
        elif reps[0] <= source_end <= reps[1]:
            len1 = len(alignment_seq[
                reps[0]:source_end].replace('-', ''))
            if len1/source_site > 0.5:
                repeat_found = True
                break

        elif source_start <= reps[0] <= source_end:
            len1 = len(alignment_seq[
                reps[0]:reps[1]].replace('-', ''))
            if len1/source_site > 0.5:
                repeat_found = True
                break
        if reps[0] <= target_start <= reps[1]:
            len1 = len(alignment_seq[
                target_start:reps[1]].replace('-', ''))
            if len1/target_site > 0.5:
                repeat_found = True
                break
        elif reps[0] <= target_end <= reps[1]:
            len1 = len(alignment_seq[
                reps[0]:target_end].replace('-', ''))
            if len1/target_site > 0.5:
                repeat_found = True
                break
        if target_start <= reps[0] <= target_end:
            len1 = len(alignment_seq[
                reps[0]:reps[1]].replace('-', ''))
            if len1/target_site > 0.5:
                repeat_found = True
                break

    return repeat_found


def get_genomic_coordinates(database,
                            collection,
                            alignment_dir,
                            seq_dir):
    """
    Get TSMs full filling the criteria from a database.
    Confirm that TSM target and source sites do not overlap with a
    repeat more than 50%. If a TSM pass the criteria, the script
    yields the information on a TSM.

    Argument
    ========

    database: A TSM databse (mongo database)
    collection: A collection name in the mongo database
    alignment_dir: A full path to a directory containing sequence alignments
    seq_dir: A directory containing non-aligned sequences.


    Yields
    ======

    human_chromosome: A chromosome
    source_start: A starting index of a TSM source site
    source_end: The last index of a TSM sources tie
    target_start: A starting index of a TSM target site
    target_end: The last index of a TSM target site
    length: the length of a TSM

    """

    from RNA_alignments import indexes_in_original_seq

    client = MongoClient()
    db = client[database]
    collection = db[collection]

    query = {"mismatches": {"$gt": 5},
             "TSM_length": {"$gt": 15},
             "quality": True}

    used_files = defaultdict(dict)
    cursor = collection.find(query, no_cursor_timeout=True, batch_size=1)
    for document in cursor:

        hit_dict = {}

        for key, value in document.items():

            hit_dict[key] = value

        TSM_source_start = int(hit_dict["tsm_source_start"])
        TSM_source_end = int(hit_dict["tsm_source_end"])
        TSM_target_start = int(hit_dict["tsm_target_start"])
        TSM_target_end = int(hit_dict["tsm_target_end"])

        human_chromosome = hit_dict["human_chromosome"]
        human_location = hit_dict["human_location"]
        identifier = hit_dict["identifier"]
        length = hit_dict["TSM_length"]
        query_node = hit_dict["query"]

        if identifier in used_files:
            if query_node in used_files[identifier]:
                if len(set(range(
                    TSM_source_start, TSM_source_end)).intersection(
                        set(range(
                            used_files[
                                identifier][query_node]["TSM_s_start"],
                            used_files[identifier][
                                query_node]["TSM_s_end"])))) >= 15:
                    continue

        bio_align_dict = SeqIO.to_dict(SeqIO.parse(
            alignment_dir + identifier + '_pagan.fas',
            "fasta"))
        align_dict = OrderedDict()

        for seq in bio_align_dict:
            if "Homo" in seq:
                Homo_sapiens = str(bio_align_dict[seq].seq)

            align_dict[seq] = str(bio_align_dict[seq].seq)
        if run_dustmasker(seq_dir + identifier,
                          Homo_sapiens,
                          TSM_source_start,
                          TSM_source_end,
                          TSM_target_start,
                          TSM_target_end) is True:
            continue

        used_files[identifier] = defaultdict(dict)
        used_files[identifier][query_node] = {"TSM_s_start": TSM_source_start,
                                              "TSM_s_end": TSM_source_end}

        aligned_indexes = [TSM_source_start, TSM_source_end,
                           TSM_target_start, TSM_target_end]

        indexes = indexes_in_original_seq(aligned_indexes,
                                          Homo_sapiens)

        chunk_start = int(human_location.split('-')[0])

        try:
            chunk_start = int(human_location.split('-')[0])
            source_start = chunk_start + indexes[0]
            source_end = chunk_start + indexes[1]
            target_start = chunk_start + indexes[2]
            target_end = chunk_start + indexes[3]
            TSM_orig_start = chunk_start + hit_dict["tsm_target_start"]
            TSM_orig_end = chunk_start + hit_dict["tsm_target_end"]
        except:
            continue

        yield (human_chromosome,
               source_start,
               source_end,
               target_start,
               target_end,
               length,
               TSM_orig_start,
               TSM_orig_end)

    cursor.close()


def write_bed_file(output,
                   chromosome,
                   start,
                   end,
                   orig_start,
                   orig_end):
    """
    Write a bedfile.


    Arguments
    ==========

    output: outputfile
    chromosome: A chromosome
    start: A starting index
    end: the last index
    """

    with open(output, 'a+') as f:
        f.write(chromosome + '\t' + str(start) + '\t' +
                str(end) + '\t' + str(orig_start) + '\t' +
                str(orig_end) + '\n')


def main():
    """
    Get TSM coordinates from a database and write human genomic
    coordinates of TSMs into two different .bed files. The other
    one contains TSM source coordinates and the other
    TSM target coordinates.
    """

    database = sys.argv[1]
    collection = sys.argv[2]
    alignment_dir = sys.argv[3]
    target_file = sys.argv[4]
    seq_dir = sys.argv[5]

    sorted_dict_target = {}

    for coordinate in get_genomic_coordinates(
            database,
            collection,
            alignment_dir,
            seq_dir):

        if coordinate[0] not in sorted_dict_target:
            sorted_dict_target[coordinate[0]] = {}

        sorted_dict_target[coordinate[0]][coordinate[3]] =\
            [coordinate[4], coordinate[6], coordinate[7]]

    order = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
             'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2',
             'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5',
             'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY']

    for chrom in order:
        sort = sorted(sorted_dict_target[chrom])
        for coord in sort:

            write_bed_file(target_file, chrom,
                           coord,
                           sorted_dict_target[chrom][coord][0],
                           sorted_dict_target[chrom][coord][1],
                           sorted_dict_target[chrom][coord][2])


if __name__ == "__main__":
    main()
