library(data.table)
library(ggplot2)
library(knitr)

setwd('/Users/nicoleferraro/Documents/Stanford/Spring_1617/STATS345/Project/graph_genome_GENE245/')

#Read in sample information file
metadata = fread('dnase_metadata_2016-12-05.txt')

#Peak overlap file with summary of all peak-variant overlaps
overlap_data = fread('dnase_data/peak_overlap.txt')
colnames(overlap_data) = c('Sample', 'Overlap', 'Total')
overlap_data$Proportion = overlap_data$Overlap / overlap_data$Total

#Variant information
variant_data = fread('data/ALL.chr2.sorted.position.vcf')
colnames(variant_data) = c('Chr', 'Start', 'End', 'RSID')
num_vars = nrow(variant_data)

#Get hypergeometric test p-value for number of overlaps seen, given below set up
#             Variants                    Non-Variants
# Peaks       overlap                     peaks without overlap       total peaks present
# Non-Peaks   variants without overlap    remainder                   total non-peak regions
#             total variants present      total non-variant regions

calc_hyper <- function(pdata) {
  overlap_peaks = as.numeric(pdata[2])
  peak_total = as.numeric(pdata[3])
  #variants + peaks
  x = overlap_peaks
  #peaks
  q = peak_total
  #non-peaks
  p = peak_total - num_vars
  #peaks
  nn = peak_total
  #P(X > x), probability that we see more than this number of peaks
  return(phyper(x, q, p, nn, lower.tail=F))
}

hyper_pvals = apply(overlap_data, 1, function(x) calc_hyper(x))
hist(hyper_pvals, main='Distribution of Hypergeomtric p-values in skin', xlab='P-value', breaks=10)

#Intersect samples with metadata information
sample_ids = gsub("_", ";", overlap_data$Sample)
sample_pvals = cbind(sample_ids, hyper_pvals)
sample_pvals_order = sample_pvals[order(sample_pvals[,1]),]
meta = metadata[which(metadata$encode_expt_accession %in% sample_ids),]
meta_order = meta[order(meta$encode_expt_accession),]

sample_data = as.data.frame(cbind(sample_pvals, meta_order$common_name, meta_order$organ, meta_order$sex, meta_order$life_stage, meta_order$description))
colnames(sample_data) = c('Sample', 'Pval', 'Common_name', 'Organ', 'Sex', 'Life_stage', 'Description')

View(sample_data)

#Look at proportion of overlaps occurring in gene region of interest
intersect_files = dir(path=paste(getwd(), '/dnase_data/', sep=''), pattern='*.intersected.bed')
abca_inters = rep(0, length(intersect_files))
totals = rep(0, length(intersect_files))
i = 1
abca_start = 214931542
abca_end = 215138591
for (inFile in intersect_files) {
  intersect_data = fread(paste('dnase_data/',inFile, sep=''))
  starts = intersect_data[,2]
  ends = intersect_data[,3]
  abca_inds = which(intersect_data[,2] >= abca_start & intersect_data[,3] <= abca_end)
  abca_inters[i] = length(abca_inds)
  totals[i] = nrow(intersect_data)
  i = i + 1
}


