#!/bin/bash

if [ ! -d "dnase_data" ]; then
  mkdir dnase_data
fi

scp nferraro@nandi.stanford.edu:/mnt/data/integrative/dnase/ENCSR*Skin*.DNase-seq/out*/peak/macs2/rep1/*pval0.1.narrowPeak.gz dnase_data/
cd dnase_data/
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
  bedtools intersect -a data/ABCA12/ALL.chr2.sorted.position.vcf -b $file.nochr.bed -wa -wb > $name.intersected.bed
done

for file in $(ls *.intersected.bed)
do
  name=$(echo $file | cut -f 1 -d '.')
  overlap=$(cat $file | wc -l)
  orig=$(cat $name.chr2.sorted.bed.nochr.bed | wc -l)
  echo $name
  echo $overlap / $orig
done
