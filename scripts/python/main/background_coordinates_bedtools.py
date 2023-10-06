"""
This script creates background coordinates.
"""
import os
import subprocess
import sys



dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def draw_coordinates(TSM_midpoint_file,
                     gene_coordinates,
                     outfile,
                     repeat):
    """
    Randomly suffles list of genomic coordinates using bedtools.
    

    Arguments
    ==========

    background_size: The background size
    length_counts: A counts of observed TSM lengths as a list format
    lengths_prob: Probabilities for different 


    Returns
    ========

    coordinates: coordinates of good quality TSMs
    """

    coordinates = {}


    p = subprocess.check_output(["bedtools",
                                 "shuffle",
                                 "-i",
                                 "{TSM_midpoint_file}".format(TSM_midpoint_file=TSM_midpoint_file),
                                 "-g",
                                 "{gene_coordinates}".format(gene_coordinates=gene_coordinates),
                                 "-incl",
                                 "{gene_coordinates}".format(gene_coordinates=gene_coordinates)])

        
    line_list = str(p, 'utf-8').split('\n')
        
   

    with open(outfile + '_' + str(repeat), 'a+') as f:

        for line in line_list:

            f.write(line + '\n')


def main():
    """
    Create background set.
    """

    outfile = sys.argv[1]
    repeat_number = sys.argv[2]
    TSM_midpoint_file = sys.argv[3]
    gene_coordinates = sys.argv[4]

    draw_coordinates(
            TSM_midpoint_file,
            gene_coordinates,
            outfile,
            repeat_number)

    
if __name__ == "__main__":
    main()
 

