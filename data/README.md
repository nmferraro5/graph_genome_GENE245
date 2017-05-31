# Data Sources

#### File with all 1000 genomes project individuals and their sex and populations assignments
1000genomes_samples_clean.csv - Reformatted
1000genomes_samples_orig.csv - Original


#### Added TAB delim format of Chr2 vcf > ALL.chr2.txt

#### Reference genome assembly used in 1000 genomes phase 3
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
samtools index human_g1k_v37.fasta.gz
```

#### Individual BAM file High Cov 1000 genomes
```
# Variables
abca12='2:215793266-216065151'
site=ftp://ftp.ncbi.nlm.nih.gov/1000genomes

# Get BAM index file
wget $site/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam.bai

# Get BAM file
wget $site/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam

# Just get region of interest (gene)
samtools view -h HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam $abca12 > HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.abca12.bam

# Use bamToCns.sh to convert BAM to Consensus file
```

#### 1000 genomes VCF files:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

```
# Get all populations subset of VCF
# Tabix is installed and vcf-subset installed on sherlock2 in /share/PI/pritch/Margaret/bin
# BGZIP only on local
tabix -h ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 17:1471000-1472000 | perl vcf-subset -c HG00098 | bgzip -c /tmp/HG00098.20100804.genotypes.vcf.gz
```
