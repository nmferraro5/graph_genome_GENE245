### PCA on 1000 genomes data

setwd("~/Documents/stanford/classes/spring17/bmi245/project/graph_genome_GENE245/")

####### Make table of samples to populations for coloring 
# VCF file just has sample IDs, no distinction for which population
# they belong to.
# Downloaded sample info list from 10000 genomes (table: igsr_samples.tsv)
# http://www.internationalgenome.org/data-portal/sample
# Removed biosampleID and dataCollection columns because missing values, not needed
info<-read.csv("data/1000genomes_samples_clean.csv")
colnames(info)[1]<-"sample.id"
key <- data.frame(sample = info$sample.id, pop = info$Population.code, superpop = info$Superpopulation.code)
data <- data.frame(id=abca12_pca$sample.id)
master <- merge(key,data,by.x = "sample", by.y = "id")
master$color <-'transparent'
master[which(master$superpop == 'AFR'), "color"]<-'black'
master[which(master$superpop == 'EUR'), "color"]<-'blue'


### VCF to PCA using SNPRelate package

source("https://bioconductor.org/biocLite.R")
biocLite("SNPRelate")
library(SNPRelate)

vcf.fn <-"data/ABCA12/ALL.chr2.214931542-215138591.phase3.genotypes.vcf"

snpgdsVCF2GDS(vcf.fn, "vcf.all",  method="biallelic.only")
genofile <- snpgdsOpen("vcf.all")
vcf_pca<-snpgdsPCA(genofile)

# By super population
plot(vcf_pca$eigenvect[,1],vcf_pca$eigenvect[,2], 
     col=as.factor(master$superpop), pch=2)
legend('topright', legend = levels(as.factor(master$superpop)), 
       xlab = "PCA - 1",
       ylab = "PCA - 2",
       col = 1:5, cex = 0.8, pch = 16)

# By population
plot(vcf_pca$eigenvect[,1],vcf_pca$eigenvect[,2], 
     col=as.factor(master$superpop), pch=2)
legend('topright', legend = levels(as.factor(master$superpop)), 
       xlab = "PCA - 1",
       ylab = "PCA - 2",
       col = 1:5, cex = 0.8, pch = 16)

# Only show EUR and AFR
plot(vcf_pca$eigenvect[,1],vcf_pca$eigenvect[,2], 
     xlab = "PCA - 1",
     ylab = "PCA - 2",
     col=master$color, pch=2)
legend('topright', legend = c("EUR","AFR"), col = c("blue","black"), 
       cex = 0.8, pch = 16)



