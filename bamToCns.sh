#!/bin/bash
# Get consensus sequence (https://www.biostars.org/p/184944/#185560)
#
#USAGE: sh bamToCns.sh <input.bam> <ref.fa>
#
#Export file paths to access samtools and bcftools and tabix
export PATH=/share/PI/pritch/Margaret/bin:$PATH
source ~/.bashrc

#TO GET REF.FA INPUT FILE
# Downlaod Fasta hg19 file used in 1k genomes 
# wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz

# TO GET INPUT.BAM INPUT FILE
# Download exome alignment file for individual in 1000k genome project
# Download script for batches on Margaret's server test directory
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/

# Get bam index file
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.bai
# Get bam file
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam
# Make sure bam and index files are in same location. 

#Access input bam file 
inputFile=$1
# Just get root name for new file naming
input=${inputFile%.*}
ref=$2

#Get variants with mpilup and bcftools call
samtools mpileup -uf $ref $input.bam | bcftools call -mv -Oz -o $input.calls.vcf.gz
tabix $input.calls.vcf.gz
cat $ref | bcftools consensus $input.calls.vcf.gz  > $input.calls.cns.fa

# Clean up

# Modify each cns.select.fa header (>....) with individual's ID and population
# Copy all sequences and paste into Kalign browser job
