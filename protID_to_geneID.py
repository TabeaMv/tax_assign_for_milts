'''
The author of this program is Tabea Merlevede.
This program takes as input a taxonomic assignment on protein level, which is the output if you use MEGAN and a
    simplified version of your GFF file.
    For advice on how to simplify the GFF file, have a look in my wiki.
The output that you recieve here is the taxonomic assignment on gene level.
That means:
- every protein ID gets subsituted with the associated gene name
- if there are at least two proteins for one gene (reason could be splicing), the protein ID in the row in which we have
    the longer protein will get substituted with the associated gene name. The rows, where there are proteins that are not
    the longest of their respective gene, get deleted.
- You recieve a file that is max. the length of your taxonomic assignment, but most likely a few rows shorter.

Your files dont change. There will be a new created file called "taxonomic_assignment_gene_level".

In the following protein_IDs have the ending "-mRNA-X" in the txt files. Protein_names are the Protein_IDs without the -mRNA-X ending.
Gene_IDs are the ones that can be found in the second column of the simplified GFF.
'''

# get the content of the simplified gff file
with open("path/to/simplified_gff",'r') as gff:
    gff_lines = gff.readlines()

gff.close()

# get the content of the taxonomic assignment on protein level
with open("path/to/taxonomic_assignment_protein_level.txt",'r+') as assign:
    assign_lines = assign.readlines()

assign.close()

# prot_to_gene contains already found longest proteins. They get substituted as soon as the outer for loop reaches the respective row.
prot_to_gene = []

# skippable_prot_IDs: already checked proteins that dont get written in the final file because the proteins were not the longest of their respective gene
skippable_prot_IDs = []

# list of final rows, later written in txt file called taxonomic assignment on gene level
final_txt_lines = []

# go through every line of the taxonomic assignment on protein level
for assign_line in assign_lines:
    #found_similar_names contains the proteins that belong to the same gene as the protein in the current row
    # found_similar_names gets resetted for every row in the taxonomic assignment
    found_similar_names = []
    # we just need the protein ID of the current row
    prot_id_assign = assign_line.split("\t")[0]

    # check if protein in the current row has already been rated as not-longest protein of a respective gene
    # if prot_ID is in skippable_prot_IDs, it gets deleted out of it and does not get written in the final txt file
    if prot_id_assign in skippable_prot_IDs:
        skippable_prot_IDs.remove(prot_id_assign)
        continue

    # check if protein in the current row has already been rated as longest protein of a respective gene
    # if prot_Id is in prot_to_gene, it gets substituted with its respective Gene_ID
    elif prot_id_assign in prot_to_gene:
        final_txt_lines.append(prot_id_assign.rsplit("-mRNA-", 1)[0] + "\t" + assign_line.split("\t")[1])
        continue

    # we want to compare the protein name, not the protein ID, so we remove the -mRNA-X ending
    prot_name_assign = prot_id_assign.rsplit("-mRNA-", 1)[0]

    # go through every entry of the gff
    for gff_line in gff_lines:
        split_gff_line = gff_line.split("\t")
        # we just need the protein ID, gene ID and the length of the protein
        prot_id_gff = split_gff_line[0]
        gene_id_gff = split_gff_line[1]
        prot_length_gff = split_gff_line[3]
        # we want to compare the protein name, not the protein ID, so we remove the -mRNA-X ending
        prot_name_gff = prot_id_gff.rsplit("-mRNA-", 1)[0]
        # find all occurences of the protein name in the gff file and add found proteins to found_similar_names
        if prot_name_assign == prot_name_gff:
            found_similar_names.append([prot_id_gff, gene_id_gff, int(prot_length_gff.split("\n")[0])])
    # get longest protein of the ones saved in found_similar_names
    longest_similar_protein = max(found_similar_names, key=lambda x: x[2])

    # remove longest protein for the respective gene from found_similar_names
    found_similar_names.remove(longest_similar_protein)

    # add every not_longest protein_ID to list that gets checked in the beginning of the outer for loop
    for x in found_similar_names:
        skippable_prot_IDs.append(x[0])

    # if the protein ID of the current row is the same as the one of the longest_similar_protein, swap it now with the
    # respective gene ID
    if prot_id_assign == longest_similar_protein[0]:
        final_txt_lines.append(longest_similar_protein[1] + "\t" + assign_line.split("\t")[1])
    # if its not the same, append it to the list that gets checked in the beginning of the outer for loop
    else:
        prot_to_gene.append(longest_similar_protein[0])

# set your desired output path here
with open('path/to/taxonomic_assignment_gene_level.txt', 'w') as taxonomic_assignment_gene_level:
    taxonomic_assignment_gene_level.writelines(final_txt_lines)

taxonomic_assignment_gene_level.close()