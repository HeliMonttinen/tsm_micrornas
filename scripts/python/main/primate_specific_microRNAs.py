"""
This script identifies the primate-specific microRNAs.
Goes through all human microRNAs and identifies the parental
sequence, which is at least 90% identical with the human microRNA.
Records if all descendants of this parental node are among primates.

In case of TSM-origin microRNAs, it calculates the sequence identity
between the query sequence and human sequence. The result is saved
into a dictionary.
"""

from collections import OrderedDict
import sys
import os

from decimal import Decimal

from Bio import SeqIO
from ete3 import Tree
from pymongo import MongoClient

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def parse_microRNAs(microRNA_file):
    """
    A microRNA file is parsed and saved
    as a dictionry.

    Arguments
    ==========

    microRNA_file: A path to a microRNA file


    Returns
    ========

    microRNA_dict: MicroRNA identifiers and locations
                   saved into a dictionary.
    """

    microRNA_dict = dict()

    with open(microRNA_file, 'r') as f:

        for line in f:

            line_splitted = line.rstrip().split()
            identifier = line_splitted[0]
            location = line_splitted[-1].rstrip()

            start = int(location.split(':')[-1].split('-')[0])
            end = int(location.split('-')[-1])
            
            if 'PATCH' in location:
                continue
    
            microRNA_dict[identifier] = dict()
            microRNA_dict[identifier][0] = start
            microRNA_dict[identifier][1] = end

    return microRNA_dict


def TSM_microRNAs(TSM_file):
    """
    Parse a file listing TSM-origin microRNAs.

    Arguments
    =========
    TSM_file: A path to TSM file

    Returns
    =======

    TSM_dictionary: TSM identifiers in a dictionary.
                    identifiers are values.
    """
    
    TSM_dict = {}

    with open(TSM_file, 'r') as f:

        for line in f:

            line_splitted = line.split('\t')
            TSM_dict[line_splitted[0]] = int(line_splitted[2])

    return TSM_dict


def screen_microRNA_loci(directory,
                         microRNA_dict,
                         TSM_dict,
                         database):
    """
    Go through microRNA loci and check if: 1) microRNA
    sequence is only among human. MicroRNA is found is 
    >seqid 90% and coverage>95%.
    """

    from RNA_alignments import indexes_in_alignment


    species_list = {"Homo_sapiens":"",
                    "Gorilla_gorilla":"",
                    "Pan_troglodytes":"",
                    "Pongo_pygmaeus":"",
                    "Nomascus_leucogenys":"",
                    "Macaca_mulatta":"",
                    "Macaca_fascicularis":"",
                    "Papio_hamadryas":"",
                    "Chlorocebus_sabaeus":"",
                    "Callithrix_jacchus":"",
                    "Saimiri_boliviensis":"",
                    "Otolemur_garnettii":""}

    primate_specific = {}
    mutations_after_TSM = {}
    for filename in os.listdir(directory):

        identifier = filename.split('/')[-1].split('.')[0]

        if not os.path.exists(directory + identifier + '_pagan.anctree'):
            continue
        tree = Tree(directory + identifier + '_pagan.anctree', format=1)
       
        bio_align_dict = SeqIO.to_dict(
                SeqIO.parse(
                    directory + identifier  + '_pagan.fas',
                    "fasta"))

        alignment_dict = OrderedDict()

        for seq in bio_align_dict:

            alignment_dict[seq] = str(bio_align_dict[seq].seq)

        for node in tree.traverse():
            if 'Homo' in node.name:
                seq_int = node.name

                try:
                    head_start = int((node.name).split('_')[-2])
                    head_end = int((node.name).rstrip().split('_')[-1])
        
                except:
                    head_end = len(alignment_dict[node.name].replace('-',''))
                    head_end = head_start + head_end

                diff_start = microRNA_dict[identifier][0] - head_start
                diff_end = len((alignment_dict[node.name]).replace('-','')) -\
                        (head_end - microRNA_dict[identifier][1])

                new_indexes = indexes_in_alignment(
                    [diff_start, diff_end], alignment_dict[node.name])

                homo_sapiens_seq = alignment_dict[seq_int][new_indexes[0]:new_indexes[1]]
                break
        
        
        if identifier not in TSM_dict:

            only_in_primates = True
    
            for seq in alignment_dict:

                if '#' in seq:
                    continue

                alignment_dict[seq] = alignment_dict[seq][new_indexes[0]:new_indexes[1]]

                ident_res = 0
                human_len = 0
                for a,b in zip(alignment_dict[seq], homo_sapiens_seq):

                    if a == b and b != '-':

                        ident_res += 1

                    if b != '-' or a != '-':

                        human_len += 1
                                    
                sequence_identity = Decimal(ident_res/human_len)
                
                if sequence_identity > 0.9:
                    species = '_'.join(seq.split('_')[:2])
                else:
                    continue
                if species not in species_list:

                    only_in_primates = False

        else:

            client = MongoClient()
            db = client[database]
            collection = db['TSMs']

            query = {"identifier":identifier,
                     "TSM_length": TSM_dict[identifier],
                     "quality": True}


            for document in collection.find(query):

                hit_dict = {}

                for key, value in document.items():
                    hit_dict[key] = value

                query = hit_dict["query"]
                only_in_primates = True
                query_seq = alignment_dict[query][new_indexes[0]:new_indexes[1]]

                ident_res = 0
                query_len = 0
                
                for a,b in zip(query_seq, homo_sapiens_seq):

                    if a == b and b != '-':

                        ident_res += 1

                    if a != '-' or b != '-':

                        query_len += 1

                sequence_identity = Decimal(ident_res/query_len)

                if sequence_identity < 1.0:

                    mutations_after_TSM[identifier] = sequence_identity
                else:
                    print(identifier)

        if only_in_primates is True:

            primate_specific[identifier] = ""

    print(len(mutations_after_TSM))
    print(len(primate_specific))


def main():
    """
    Run the scripts
    """

    microRNA_file = sys.argv[1]
    TSM_file = sys.argv[2]
    directory = sys.argv[3]
    database = sys.argv[4]

    microRNA_dict =  parse_microRNAs(microRNA_file)        

    TSM_dict = TSM_microRNAs(TSM_file)

    screen_microRNA_loci(directory,
                         microRNA_dict,
                         TSM_dict,
                         database)

if __name__ == "__main__":
    main()
