from collections import defaultdict

line_dict = defaultdict(list)

with open("results/mouse_mirbasedb2/list_mouse_mir.txt", "r") as f:

    for line in f:

        line_splitted = line.split()
        help_list = [int(line_splitted[2]) - 50000,
                     int(line_splitted[3]) + 50000]

        line_dict["chr" + line_splitted[1].rstrip()].append(help_list)


with open(
    "data/Mus_musculus.GRCm39.109.gff3", "r"
) as f2:

    found = False
    test = False
    name_dict = dict()
    for line in f2:

        line_splitted = line.split("\t")
        if len(line_splitted) < 7:
            continue

        if "gene" in line_splitted[2]:

            if "Name" in line:
                Name = line_splitted[8].split(";Name=")[-1].split(";")[0]
            else:

                bio = line_splitted[8].split("biotype=")[-1].split(";")[0]
                gene_id = line_splitted[8].split("ID=gene:")[-1].split(";")[0]

                Name = bio + "_" + gene_id

            if "chr" + line_splitted[0].rstrip() in line_dict:

                for a in line_dict["chr" + line_splitted[0]]:

                    if a[0] <= int(line_splitted[3]) <= a[1]:

                        found = True
                        print("###")
                        break

                    elif a[0] <= int(line_splitted[4]) <= a[1]:

                        found = True
                        print("###")
                        break

                    elif int(line_splitted[3]) <= a[0]\
                            <= int(line_splitted[4]):
                        found = True
                        print("###")
                        break

                    elif int(line_splitted[3]) <= a[1] <=\
                            int(line_splitted[4]):
                        found = True
                        print("###")
                        break

                    else:
                        found = False
                        test = False

        if "biological" in line_splitted[2]:
            found = False

        if found is True:

            print("chr" + line.rstrip())

            if "ID=transcript:" in line:
                trans = line_splitted[8].split(
                        "ID=transcript:")[1].split(";")[0]
                name_dict[trans] = Name

print("###")


with open("results/mouse_mirbasedb2/transcript_symbols_mouse.txt", "w") as f:

    keys = sorted(name_dict.keys())

    for key in keys:

        f.write(key + "\t" + name_dict[key] + "\n")
