"""
This script downloads the alignments using
a given list of Ensebl identifiers and locations.
"""
import os
import sys

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def parse_Ensembl_file(data):
    """
    Yield the location of miRNA (and -+100 nt)
    in the sequence.
    """
    used = set()
    with open(data, 'r') as f:

        for line in f:

            splitted_line = line.rstrip().split()
            if splitted_line[0] in used:
                continue
            chromosome =  splitted_line[0]
            start = int(splitted_line[1])-50
            end = int(splitted_line[2])+50
            accession = splitted_line[3]

            yield (accession, chromosome + ':' + str(start) + '-' + str(end))


def main():
    """
    Download the alignments based on the
    Ensembl identifier and sequence location.
    """

    from ensembl_download import request_tool

    data = sys.argv[1]
    species = sys.argv[2]
    ensembl_alignment = sys.argv[3]
    target_dir_alignment = sys.argv[4]
    target_dir_tree = sys.argv[5]

    for accession, location in parse_Ensembl_file(data):

        try:
            obj = request_tool(species, location, ensembl_alignment)
            with open(target_dir_tree+ accession + '.tree', 'w') as f:
                f.write(obj["tree"])
            with open(target_dir_alignment + accession, 'w') as f:
                for item in obj["alignments"]:
                    if item["strand"] == 1:
                        strand = '_+_'
                    else:
                        strand = '_-_'
                    if 'Ancestor' in  item["seq_region"]:
                        continue
                    f.write(">" + item["species"] + '_' + str(item["seq_region"]) + '_' + str(item["start"]) +
                            '_' + str(item["end"]) + strand + '\n')
                    f.write(item["seq"] + '\n')
        except:
            print(accession, "failed")


if __name__ == "__main__":
    main()
