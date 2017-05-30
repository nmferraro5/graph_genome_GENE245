#!/bin/bash

if [ ! -d "dnase_data" ]; then
  mkdir dnase_data
fi

python add_pos_vcf.py 'data/ALL.noheader.chr2.txt' 'data/ALL.chr2.position.vcf'
bedtools sort -i data/ALL.chr2.position.vcf > data/ALL.chr2.sorted.position.vcf
rm data/ALL.chr2.position.vcf

cd dnase_data
scp nferraro@nandi.stanford.edu:/mnt/data/integrative/dnase/ENCSR*Skin*.DNase-seq/out*/peak/macs2/rep1/*pval0.1.narrowPeak.gz .
gunzip *.gz

for file in $(ls *.narrowPeak)
do
  name=$(echo $file | cut -f 1 -d '.')
  grep "chr2" $file | bedtools sort > $name.chr2.sorted.bed
  rm $file
done

for file in $(ls *.bed)
do
  name=$(echo $file | cut -f 1 -d '.')
  cat $file | sed 's/^chr//' > $file.nochr.bed
  rm $file
  bedtools intersect -a ../data/ALL.chr2.sorted.position.vcf -b $file.nochr.bed -wa -wb > $name.intersected.bed
done

output_file='peak_overlap.txt'
for file in $(ls *.intersected.bed)
do
  name=$(echo $file | cut -f 1 -d '.')
  overlap=$(awk -F '\t' '{print $8}' $file | sort | uniq -c | wc -l)
  orig=$(cat $name.chr2.sorted.bed.nochr.bed | wc -l)
  echo -e $name"\t"$overlap"\t"$orig >> $output_file
done
