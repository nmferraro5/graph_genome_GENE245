### Margaret Antonio
### GENE245 Project with Nicole Ferraro

### Functions to run PCA on 1000 genomes data for two populations
### Requires getPops() function in knnOnGenotypes.R
### Used for Fig.1 Panels C and D

# Packages required
library(ggfortify)

# Data file paths
# These files were downloaded from 1000 genomes site and formatted 
# using get_tkg_data.sh
home="~/Documents/stanford/classes/spring17/bmi245/project/graph_genome_GENE245/data/"
abca12GenoFile = paste0(home,"ABCA12/ALL.2.215796266-216003151.justGenotypes.txt")
znfGenoFile = paste0(home,"ZNF/ALL.2.185463093-185804214.justGenotypes.txt")
randomGenoFile = paste0(home,"CHR2/ALL.chr2.biallelic.genotypes.random6k.txt")

# Data frames for each set
tgp.znf <- getPops(znfGenoFile,'EAS','AFR')
tgp.rand <- getPops(randomGenoFile,'EUR','AFR')
tgp.abca12 <- getPops(abca12GenoFile,'EUR','AFR')

# Plot PCA
plotPCA <- function(data, title){
  df.data <- data[c(1:ncol(data)-1)]
  autoplot(prcomp(df.data),data = data, colour = 'superpop', main = title)
}

### Figure 1, Panels C and D

# Panel C
plotPCA(tgp.rand, "African and European Populations; 6000 Random Chromosome 2 variants")

# Panel D
plotPCA(tgp.abca12, "African and European Populations; ABCA12 variants")

