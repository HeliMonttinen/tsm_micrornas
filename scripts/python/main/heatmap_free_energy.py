"""
Create a presence-absence matrix colored based on
the free-energy changes. Only those microRNAs are
included, which can be fully explained by TSM.
"""

import os
import sys
from collections import defaultdict
from decimal import Decimal
from ete3 import Tree
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.font_manager import FontProperties
import numpy as np
import pandas as pd
import seaborn as sns


dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def microRNA_int(microRNA_file):
    """
    Parses a microRNA file containing microRNA
    identifiers and information. Returns a dictionary
    where Ensembl identifier is a key and a value
    is a microRNA identifier. If multiple microRNAs are
    processed from the same locus, an asterisk is
    added to the end.

    Argument
    ========
    microRNA_file: A path to the microRNA info file

    Output
    =======
    microRNAs: A dictionary, where microRNA is an
               identifier and microRNA name as a value
    """

    microRNAs = dict()
    with open(microRNA_file) as f:

        for line in f:

            line_splitted = line.split('\t')
            if len(line_splitted[4].split(',')) > 1:
                micronames = '-'.join(
                        line_splitted[5].lstrip().split(
                            ',')[0].rstrip('a').rstrip(
                                'b').split('-')[:3]) + '*'

            else:
                micronames = '-'.join(
                        line_splitted[5].strip().rstrip(
                            'a').rstrip('b').split('-')[:3])

            microRNAs[line_splitted[0]] = micronames

    return microRNAs


def collect_information_taxa(microrna_int_list,
                             alignment_dir,
                             tree_file):
    """
    Collect information on the species in a given
    phylogenetic tree. Includes only those species,
    which has are present at least in 10 alignments.
    Removes those taxa, which are not included the species list.

    Arguments
    =========
    microrna_inf_list: An ouput from microRNA_inf function
    alignment_dir: A path to the alignment directory
    tree_file: The full-size phylogenetic species containing
               all species

    Returns
    =======
    species_list: Species to be included in the heatmap
    tree: Phylogenetic tree containing only the species
          in the species list
    """

    species_set = dict()
    for microRNA in microrna_int_list:
        pathname = alignment_dir + microRNA + '.fas'

        with open(pathname, 'r') as f:

            for line in f:
                if '>' in line and '#' not in line:

                    tax_name = line.lstrip('>').rstrip().split('_')
                    count = 0
                    for word in tax_name:
                        number = False

                        for char in word:
                            if char.isnumeric():
                                number = True
                        if number is True:
                            break
                        elif not (word.startswith('chr') or
                                  word.startswith('scaffold')):
                            count += 1

                    if tax_name[0].islower():

                        if tax_name[1].islower():
                            species_name = ' '.join(tax_name[2:count])
                        else:
                            species_name = ' '.join(tax_name[1:count])
                    else:
                        species_name = ' '.join(tax_name[:count])

                    if species_name not in species_set:
                        species_set[species_name] = 1
                    else:
                        species_set[species_name] += 1

    leaf_nodes_present = set()
    leaf_nodes_not_present = set()
    species_list = []
    tree = Tree(tree_file)

    leaf_nodes = tree.get_leaves()

    for leaf in leaf_nodes:
        found = False
        leaf_name = leaf.name
        spes = leaf_name.split('_')[:2]
        full_name = ' '.join(spes)

        for spes1 in species_set:
            if full_name in spes1:
                leaf_nodes_present.add(leaf_name)
                species_list.append(spes1)
                found = True

        if found is False:
            leaf_nodes_not_present.add(leaf_name)

    nodes_to_detach = set()
    for node in tree.traverse():
        if node.is_leaf():
            continue
        else:
            leafes = node.get_leaves()

            present = False
            for sub_leaf in leafes:

                if sub_leaf.name in leaf_nodes_present:

                    present = True
            if present is False:
                nodes_to_detach.add(node)

    for node in nodes_to_detach:
        node.detach()

    for node in tree.get_leaves():
        if node.name not in leaf_nodes_present:
            node.detach()
        else:
            node.name = (node.name).replace('_', ' ')

    return species_list, tree


def parse_parse_temperatures(structure_dir, species_list):
    """
    Parses minimum free energies from .db (dot-parenthesis)
    for all species in each loci.

    Arguments
    =========

    structure_dir: A directory containing strutcure files
    species_list: A list of species to be included in the heatmap

    Returns
    =======

    free_energy_dict: Free energies as a dictionary

    """

    free_energy_dict = defaultdict(dict)

    for subdir, dirs, files in os.walk(structure_dir):

        for file in files:
            if '#' in file:
                continue
            non_struct = False
            if not file.endswith('.db'):
                filename = file.rstrip('.fasta')
                filename = filename + '.db'
                filepath = structure_dir + filename
                if not os.path.exists(structure_dir + filename):
                    non_struct = True
            else:
                continue

            microRNA = file.split('_')[0]

            tax_name = file.rstrip().split('_')[1:]
            count = 0
            for word in tax_name:
                number = False

                for char in word:
                    if char.isnumeric():
                        number = True
                if number is True:
                    break
                elif not (word.startswith('chr') or
                          word.startswith('scaffold')):
                    count += 1

            if tax_name[0][0].islower():

                if tax_name[1][0].islower():
                    species = ' '.join(tax_name[2:count])
                else:
                    species = ' '.join(tax_name[1:count])
            else:
                species = ' '.join(tax_name[:count])

            found = False
            for spes in species_list:
                if species in spes:

                    species_name = spes
                    found = True
            if found is False:
                continue

            if non_struct is False:
                with open(filepath, 'r') as f:

                    try:
                        last_line = f.readlines()[-1]
                        free_energy_dict[microRNA][species_name] = 0
                    except:
                        continue

                    free_energy = Decimal(
                            last_line.rstrip().rstrip(')').split(' (')[-1])

                free_energy_dict[microRNA][species_name] = free_energy
            else:

                free_energy_dict[microRNA][species_name] = 0
    return free_energy_dict


def organise_columns(TSM_table, tree, species_list, microRNAs_int):
    """
    Organises microRNA identifiers into columns based on
    the phylogenetic depth of the TSM events.

    Arguments
    ==========

    TSM_table: A TSM table containing information about
               TSMs and which species the loci are found.
    tree: A phylogenetic tree file containing all species used
          for a heatmap.
    species_list: A list of species to be included in the heatmap
    microRNAs_int: The microRNAs of interest

    """

    column_order = []

    microRNA_species = defaultdict(list)
    size_list = set()
    size_list_dict = defaultdict(list)

    with open(TSM_table, 'r') as f:

        for line in f:

            line_splitted = line.rstrip().split('\t')

            microRNA_id = line_splitted[0]
            species = line.split('\t')[20]
            species_in_table = species.split(', ')
            species_present = []
            species_present_nodes = []

            for spes1 in species_in_table:
                spes1 = spes1.rstrip().replace('_', ' ')

                for spes2 in species_list:

                    if spes2[0].islower():
                        if not spes2[1].islower():
                            spes2 = ' '.join(spes2.split()[1:])
                        else:
                            spes2 = ' '.join(spes2.split()[2:])

                    if spes1 in spes2:

                        species_present.append(spes2)
                        species_present_nodes.append(tree&spes2)

            if len(species_present_nodes) > 1:
                common_ancestor = tree.get_common_ancestor(
                        *species_present_nodes)

                for leaf in common_ancestor.iter_leaves():

                    species_present.append(leaf.name)
            else:
                species_present.append(species_present_nodes[0].name)

            species_len = len(list(set(species_present)))

            if microRNA_id in size_list_dict and\
                    species_len > size_list_dict[microRNA_id][0]:
                microRNA_species[
                    size_list_dict[microRNA_id][0]].remove(microRNA_id)
                microRNA_species[species_len].append(microRNA_id)

            elif microRNA_id not in size_list_dict:
                microRNA_species[species_len].append(microRNA_id)

            size_list.add(species_len)
            size_list_dict[microRNA_id].append(species_len)
    size_list = sorted(list(size_list), reverse=True)
    column_names = []
    for size in size_list:
        hsa_ids = []
        nums = []
        for microRNA_id in microRNA_species[size]:

            if microRNA_id in microRNAs_int:
                hsa_ids.append(microRNAs_int[microRNA_id])
                num = int(microRNAs_int[microRNA_id].split('-')[2].rstrip('*'))
                nums.append(num)
        nums = sorted(nums)
        for num in nums:
            for hsa in hsa_ids:
                hsa_num = hsa.split('-')[2].rstrip('*')
                if str(num) == hsa_num:
                    column_names.append(hsa.lstrip())
                    for a in microRNAs_int:
                        if microRNAs_int[a] == hsa:
                            column_order.append(a)

    return column_order, size_list_dict, column_names


def create_table(
        species_list,
        free_energy_dict,
        column_order,
        outputfile,
        size_len_dict,
        column_names):
    """
    A script to create a heatmap.

    Arguments
    =========

    species_list: Species names to be included the heatmap
    free_energy_dict: A nested dictionary, where a microRNA
                      identifier is the first key, species name is
                      the second key and a free energy is a value.
    column_order: The list, where microRNA identifiers indicate
                  the order for the columns.
    outputfile: The outputfile where
    size_len_dict: Dictionary indicating how many species are
                    affected by a TSM
    column_names: Column names (a list)

    """
    rows = species_list

    columns = defaultdict(list)
    minimum = 0
    maximum = -1000

    new_free_energy_dict = defaultdict(dict)
    for ids in free_energy_dict:
        minimum = -1000
        maximum = 0
        for spes in free_energy_dict[ids]:
            if free_energy_dict[ids][spes] < maximum:
                maximum = free_energy_dict[ids][spes]
            elif free_energy_dict[ids][spes] > minimum:
                minimum = free_energy_dict[ids][spes]
        for spes in free_energy_dict[ids]:
            new_free_energy_dict[ids][spes] =\
                    free_energy_dict[ids][spes]/(-maximum)

    for spes in rows:

        for col in column_order:

            if col in new_free_energy_dict:

                if spes in new_free_energy_dict[col]:

                    if new_free_energy_dict[col][spes] > maximum:
                        maximum = free_energy_dict[col][spes]
                    elif new_free_energy_dict[col][spes] < minimum:
                        minimum = free_energy_dict[col][spes]

                    columns[col].append(
                        round(float(new_free_energy_dict[col][spes]), 2))
                else:
                    columns[col].append(np.nan)
            else:
                columns[col].append(np.nan)

    df = pd.DataFrame(columns, columns=column_order, index=species_list)

    df.to_csv(path_or_buf="heat_temperature.csv", sep='\t')

    f, ax = plt.subplots(figsize=(30, 15))

    cmap = sns.color_palette("viridis", as_cmap=True)
    sns.heatmap(df.isnull(), cmap=['lightgrey'], cbar=False)

    ax = sns.heatmap(df, vmin=-1, vmax=0, square=False,
                     cbar=True, cmap=cmap, center=-0.5)
    plt.gcf().axes[1].invert_yaxis()

    for i in range(len(column_order)):

        for a in range(len(size_len_dict[column_order[i]])):
            minimum = min(set(size_len_dict[column_order[i]]))
            if len(set(size_len_dict[column_order[i]])) == 1:
                rect = patches.Rectangle((i, 0),
                                         1.0,
                                         size_len_dict[column_order[i]][a],
                                         linewidth=2.0,
                                         edgecolor='black',
                                         facecolor='none')
            elif size_len_dict[column_order[i]][a] == minimum:
                rect = patches.Rectangle((i+0.1, 0),
                                         0.8,
                                         size_len_dict[column_order[i]][a],
                                         linewidth=2.2,
                                         edgecolor='red',
                                         facecolor='none')
            else:
                rect = patches.Rectangle((i, 0),
                                         1.0,
                                         size_len_dict[column_order[i]][a],
                                         linewidth=2.0,
                                         edgecolor='black',
                                         facecolor='none')

            ax.add_patch(rect)

    plt.tick_params(axis='x',
                    which='major',
                    labelrotation=90,
                    labelbottom=False,
                    bottom=False,
                    top=False,
                    labeltop=True)

    species_font = FontProperties()
    species_font.set_name('Arial')
    species_font.set_style('italic')
    species_font.set_size('medium')
    ax.set_xticklabels(column_names, size="large")
    ax.set_yticklabels(rows, fontstyle="italic", size="large")
    f.set_tight_layout(True)

    plt.savefig(outputfile)


def main():
    """
    Run functions for creating a heatmap.
    """

    microRNAs_table = sys.argv[1]
    tree_file = sys.argv[2]
    alignment_dir = sys.argv[3]
    structure_dir = sys.argv[4]
    outputfile = sys.argv[5]
    output_tree = sys.argv[6]

    microRNAs_int = microRNA_int(microRNAs_table)

    species, new_tree = collect_information_taxa(
            microRNAs_int,
            alignment_dir,
            tree_file)

    free_energy_dict = parse_parse_temperatures(
            structure_dir,
            species)

    column_order, size_list_dict, col_names = organise_columns(
            microRNAs_table,
            new_tree,
            species,
            microRNAs_int)

    create_table(
        species,
        free_energy_dict,
        column_order,
        outputfile,
        size_list_dict,
        col_names)

    new_tree.write(outfile=output_tree)


if __name__ == "__main__":
    main()
