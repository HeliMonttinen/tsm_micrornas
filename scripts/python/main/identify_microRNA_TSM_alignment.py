"""
"""

import sys



def identify_mir_locus(mirs):
    """
    """
    identifiers = {}
    bad_count = 0

    with open(data_file, 'r') as f:
        for line in f:
            line_splitted = line.rstrip().split()
            mir_identifier = line_splitted[3]
            if not mir_identifier.endswith('_pri'):
                continue

            chromosome = line_splitted[0]

            start = int(line_splitted[1])
            end = int(line_splitted[2])
            reg = set(range(start, end))
            found_same_seq = False

            for loc in locus_dict:
                chr_set = loc.split(':')[0]
                if chr_set == chromosome:
                    start_set = int(loc.split(':')[1].split('-')[0])
                    end_set = int(loc.split(':')[1].split('-')[1])
                    set_reg = set(range(start_set, end_set))
                    intersec = reg.intersection(set_reg)
                    if len(intersec) > 0:
                        print("bad", mir_identifier)
                        found_same_seq = True
                        bad_count += 1
                        break

            if found_same_seq is False:
                locus = chromosome + ':' + str(start) + '-' + str(end)
                locus_dict[locus] = mir_identifier

    for locs in locus_dict:
        chr_n = str(locs.split(':')[0])
        start = int((locs.split(':')[-1]).split('-')[0])
        end = int(locs.split('-')[-1])
        identifiers[locus_dict[locs]] = [chr_n, start, end]

    print(bad_count)
    return identifiers


def main():
