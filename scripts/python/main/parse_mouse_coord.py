from liftover import get_lifter
import sys
from collections import defaultdict


mouse_file = sys.argv[1]

microRNA_name = 'None'
products = {}
micro_dict = defaultdict(dict)

line_number = 0
line_number_dict = defaultdict(list)


converter = get_lifter('mm10', 'mm39')

with open(mouse_file, 'r') as f:
    line_number = 0
    for line in f:
        
        if line.startswith('chr') and 'primary' in line:
            line_splitted = line.rstrip().split('\t')
            chrom = line_splitted[0].lstrip('chr')
            print(chrom, line_splitted[3].rstrip())
            
            start = converter[chrom][int(line_splitted[3])][0][1]
            end = converter[chrom][int(line_splitted[4])][0][1]

            strand = line_splitted[6]

            other_info = line_splitted[-1].split(';')

            for info in other_info:

                if 'mmu' in info:

                    microRNA_name = '-'.join(info.split('Name=')[-1].rstrip().lower().split('-')[:3])
                    orig = info.split('Name=')[-1]
                    products[microRNA_name] = dict()

                if 'ID' in info:

                    microname2 = info.split('=')[-1]

            micro_dict[microRNA_name] = {"chromosome": chrom,
                                         "start": start,
                                         "end": end,
                                         "strand": strand,
                                         "microname2": microname2}
        low = microRNA_name.lower()
        if low in line.lower() and 'MIMAT' in line:
            
            info_mimat = line.split('\t')[-1].split(';')
            prod_chrom = line_splitted[0].lstrip('chr')
            line_splitted = line.rstrip().split('\t')
            prod_start = converter[prod_chrom][int(line_splitted[3])][0][1]
            prod_end = converter[prod_chrom][int(line_splitted[4])][0][1]

            for info2 in info_mimat:
                if 'ID' in info2:

                    product_id = info2.split('=')[-1]

                    products[low][product_id] = {'start': prod_start,
                                                 'end': prod_end}

with open('results/mouse_mirbasedb2/mouse_mir_coords.txt', 'w') as f2:
    with open('results/mouse_mirbasedb2/microRNA_mir_info.txt', 'w') as f3:

        for item in micro_dict:

            f2.write(micro_dict[item]['chromosome'] + '\t' + str(micro_dict[item]['start']) + '\t' + str(micro_dict[item]['end']) + '\t' + item +'\n')

            f3.write(micro_dict[item]['chromosome'] + ':' + str(micro_dict[item]['start']) +'-' + str(micro_dict[item]['end']) + '\t' + item +
                    '\t' + micro_dict[item]["microname2"] + '\t' + micro_dict[item]["strand"] + '\t' )
            
            prod_line = ""
            text = ""
            
            for prod in products[item]:

                text = text + prod+'('+micro_dict[item]['chromosome'] +':' + str(products[item][prod]["start"]) + '-' + str(products[item][prod]["end"]) + ');'

            text=text.rstrip(';')
            f3.write(text+ '\n') 
