# balancing selection window popgenome
bedtools intersect -a Data/Divergent_Masks/Masks_By_Type/Window_Balancing_selection.bed -b Data/Ce_Annotations/gene.protein_coding.bed -wo |\
cut -f -8,14 |\
cut -f 1,4 -d';' |\
cut -f 1,3 -d'=' |\
sed 's/ID=/CELE_/g'|\
sed 's/Gene://g' |\
awk -F"\t" '!seen[$5, $6, $7, $8]++' > Data/Divergent_Masks/Masks_By_Type/Window_Balancing_Selection_Genes.bed