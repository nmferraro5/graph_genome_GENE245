# Data for Population Variants for ABCA12


#### 1. Download phase 3 all variants for chr2
Phase 3 (26 populations, 2504 individuals) CHR2 
Download link-> ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
Download site-> ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

#### 2. Get ABCA12 variants only
`tabix -h ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 2:214931542-215138591 | bgzip -c > ALL.chr2.214931542-215138591.phase3.genotypes.vcf.gz`

Next time, pipe it
`tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 2:214931542-215138591 | bgzip -c > ALL.chr2.214931542-215138591.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`

`tabix -h link_to_input_vcf.gz chr:start-end | bgzip -c > input_region_vcf.gz`

#### 3. Check stats for vcf

`rtg vcfstats vcfinput.vcf.gz | bgzip -c > vcfinput.stats.txt`
ALL.chr2.214931542-215138591.phase3.stats.txt
Nothing really interseting because all by individual, not super population

#### 4. Just get the SNP info not the 0|0 genotype for every single individual. Just want AFs
grep region upstream and downstream AF info field
Resulting file: 
ALL.chr2.214931542-215138591.phase3.genotypes.justPop.clean.txt
A little messed up but most SNPs are there

#### 5. Get number of multi-allelic SNPs in region 
Multi-allelic SNPs have "," separation...only commas in VCF

grep , vcfinput.vcf | wc -c 
30 multiallelic SNPs in CHR2 ABCA12



## 10000 genomes project ALLELE FREQ INFO

Source: 1000 genomes README for call variants (ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_vcf_info_annotation.20141104)

The 2504 samples in the phase3 release are from 26 populations which can be categorised into five super-populations 
by continent (listed below).  As well as the global AF in the INFO field. We added AF for each super-population to the INFO field.

East Asian	EAS
South Asian	SAS
African		AFR
European	EUR
American	AMR

These allele frequences were calculated by counting the AC and AN for all the individuals from a particular super population and using that
to calculate the AF. The info tag which represents the AFs are EAS_AF, EUR_AF, AFR_AF, AMR_AF and SAS_AF

The super population assignment for each sample can be found in integrated_call_samples_v3.20130502.ALL.panel

AF for multi-allelic variants are reported for each allele independently, separated by ",".
