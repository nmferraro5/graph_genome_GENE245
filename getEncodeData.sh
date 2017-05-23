#!/bin/bash

#Files pre-chosen from online browser to select samples with appropriate data

#DNase-seq Homo sapiens uterus female adult (51 year)
wget https://www.encodeproject.org/files/ENCFF645GPP/@@download/ENCFF645GPP.bed.gz

#DNase-seq Homo sapiens transverse colon male adult (37 years)
wget https://www.encodeproject.org/files/ENCFF971UGR/@@download/ENCFF971UGR.bed.gz

#DNase-seq Homo sapiens B cell female adult (43 years)
wget https://www.encodeproject.org/files/ENCFF001WEP/@@download/ENCFF001WEP.bed.gz

gunzip *.bed.gz

#Subset files to chr7
for file in $(ls *.bed)
do
  name=$(echo $file | cut -f 1 -d '.')
  grep "chr7" $file > $name.chr7.bed 
  rm $file
done

#Insert bedtools multiple intersect all files?
