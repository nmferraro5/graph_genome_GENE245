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
