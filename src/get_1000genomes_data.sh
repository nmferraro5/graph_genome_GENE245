# Data for Population Variants for ABCA12

#### 1. Download phase 3 all variants for chr2
#Phase 3 (26 populations, 2504 individuals) CHR2

GENO_LINK="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
CHR="2"
#REGION="215796266-216003151" #ABCA12
REGION="185463093-185804214" #ZNF804A

tabix -h $GENO_LINK $CHR:$REGION | \
bgzip -c > ALL.$CHR.$REGION.genotypes.vcf.gz

#Extract Individuals and RSIDS. Convert to TAB DELIM.
# Change genotype notation to numeric

gzcat ALL.$CHR.$REGION.genotypes.vcf.gz | bcftools view --max-alleles 2 |\
sed '/##/d' | sed 's/;/        /g' | \
cut -f3,10- | sed -e 's/0|0/0/g' | sed -e 's/1|0/1/g' | sed -e 's/0|1/1/g' | \
sed -e 's/1|1/2/g' > ALL.$CHR.$REGION.justGenotypes.txt



###### 2. Download full chromosome, perform same tab delim conversion
# then randomly sample same number of variants

tabix -h $GENO_LINK $CHR | bcftools view --max-alleles 2 | sed '/##/d' | (head -n 2   && tail -n +3  | gshuf -n 6000) sed 's/;/        /g' | cut -f3,10- | sed -e 's/0|0/0/g' | sed -e 's/1|0/1/g' | sed -e 's/0|1/1/g' | sed -e 's/1|1/2/g' |   \
    (head -n 2   && tail -n +3  | gshuf -n 6000) > ALL.$CHR.random.txt


####### 3. Download population assignment of individuals
# Have to visit site and click download link http://www.internationalgenome.org/data-portal/sample
cut -f1,6 igsr_samples.tsv > 1000genomes_samples_pop.tsv
