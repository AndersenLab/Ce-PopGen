#!/bin/bash

# common 
bedtools intersect -a Data/Divergent_Masks/Masks_By_Type/Common_Freq_Masks_All_Classes.bed -b Data/Ce_Annotations/gene.protein_coding.bed -wo |\
cut -f -8,14 |\
cut -f 1,4 -d';' |\
cut -f 1,3 -d'=' |\
sed 's/ID=/CELE_/g'|\
sed 's/Gene://g' |\
awk -F"\t" '!seen[$5, $6, $7, $8]++' > Data/Divergent_Masks/Masks_By_Type/Common_Mask_Genes.bed

# intermediate 
bedtools intersect -a Data/Divergent_Masks/Masks_By_Type/Intermediate_Freq_Masks_All_Classes.bed -b Data/Ce_Annotations/gene.protein_coding.bed -wo |\
cut -f -8,14 |\
cut -f 1,4 -d';' |\
cut -f 1,3 -d'=' |\
sed 's/ID=/CELE_/g'|\
sed 's/Gene://g' |\
awk -F"\t" '!seen[$5, $6, $7, $8]++' > Data/Divergent_Masks/Masks_By_Type/Intermediate_Mask_Genes.bed

# rare 
bedtools intersect -a Data/Divergent_Masks/Masks_By_Type/Low_Freq_Masks_All_Classes.bed -b Data/Ce_Annotations/gene.protein_coding.bed -wo |\
cut -f -8,14 |\
cut -f 1,4 -d';' |\
cut -f 1,3 -d'=' |\
sed 's/ID=/CELE_/g'|\
sed 's/Gene://g' |\
awk -F"\t" '!seen[$5, $6, $7, $8]++' > Data/Divergent_Masks/Masks_By_Type/Low_Mask_Genes.bed
