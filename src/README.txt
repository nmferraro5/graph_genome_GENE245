Name: get_dnase_data.sh
Author: Nicole Ferraro
Purpose: Download all peak bed files for a given tissue from Nandi, re-format data and intersect with variants
Inputs: Name of tissue (i.e. Skin, Brain), Cluster username
Outputs: Subset chr2 bed files in dnase_data directory, *.intersected.bed files for each peak file intersected with 1000 Genomes vcf, overlap_peaks file with summary of intersections
Sample call: ./get_dnase_data.sh Skin nferraro
Requires: bedtools installed locally and on path, Nandi access

Name: add_pos_vcf.py
Author: Nicole Ferraro
Purpose: Reformat vcf file to have start and end position for using bedtools intersect
Inputs: VCF file at specified location
Outputs: New vcf file, subset to only chromosome, start, end, and rsID
Sample call: python add_pos_vcf.py
Requires: VCF file

Name: peak_stats.R
Author: Nicole Ferraro
Purpose: Script to analyze peak data, and variant intersection data, produces plots of variant-peak distributions in both tissues
Inputs: Files created from get_dnase_data.sh in the dnase_data folder, and the 1000 Genomes vcf file
Outputs: Several statistical test outcomes described in write-up, histogram showing distribution of peaks and variants in different regions (not automatically saved)
Sample call: Rscript peak_stats.R
Requires: data.table, fitdistrplus, and ggplot2 libraries, all data files produced above


Name: get_1000genomes_data.sh
Author: Margaret Antonio
Purpose: Script to download and reformat 1000 genomes data for given regions and randomly subsampling a chromosome
Inputs: None. Unless other chromosome/gene region or size of random sample is to be specified
Outputs: Tab delimited genotype file with genotypes numerically encoded
Sample call: sh get_1000genomes_data.sh

Name: knnOnGenotypes.R         
Author: Margaret Antonio
Purpose performs KNN on genotype data. Includes functions to
              (a) Read genotype file as data frame for two classes
              (b) Train KNN model 
              (c) Output ROC curves
All functions to create figures of the paper (panels A and B in Fig. 1) are called in the script
Inputs: Genotype files from get_1000genomes_data.sh. Specify file paths in the function calls at the bottom of the script.
Outputs: Two ROC plots, one for random variants in chromosome 2 and another for random variants in ABCCA12
Sample call: Rscript knnOnGenotypes.R
Requires: data files produced by get_1000genomes_data.sh, caret, e1071, pROC

Name: pcaOnGenotypes.R 
Author: Margaret Antonio
Purpose: runs PCA on genotype files
Inputs: genotype files created by get_1000genomes_data.sh and read in by function getPops() in knnOnGenotypes.R
Output: PCA plot for populations (Fig.1 Panels C and D)
Sample call: Rscript pcaOnGenotypes.R
Requires: getPops() function in knnOnGenotypes.R, ggfortify
