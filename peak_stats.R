library(data.table)
library(fitdistrplus)
library(ggplot2)

setwd('/Users/nicoleferraro/Documents/Stanford/Spring_1617/STATS345/Project/graph_genome_GENE245/')

###Peak enrichment analysis

#Function to calculate peak enrichment in a gene region as compared to background
#Inputs: list of sample peak files, gene start and end location, and a boolean if plots are created
#Returns a list of p-vals describing significance of peak enrichment for that region and peak counts
get_peak_enrichment <- function(peak_files, gene_start, gene_end, gene_plot)
{
  gene_peaks = rep(0, length(peak_files))
  gene_pvals = rep(0, length(peak_files))
  totals = rep(0, length(peak_files))
  gene_length = gene_end - gene_start
  i = 1
  #Consider parallelizing if number of samples increases drastically in future
  for (inFile in peak_files) {
    peak_data = fread(paste('dnase_data/',inFile, sep=''))
    starts = peak_data[,2]
    ends = peak_data[,3]
    gene_inds = which(peak_data[,2] >= gene_start & peak_data[,3] <= gene_end)
    #generate 1000 random regions of same length to sample background distribution of peaks
    rand_inds = sample(min(starts):max(ends)-gene_length, 1000) 
    bg_peaks = unlist(lapply(rand_inds, function(x) length(which(peak_data[,2] >= x & peak_data[,3] <= x+gene_length))))
    #Fit a negative binomial
    exp_fit = fitdist(bg_peaks, "nbinom")
    hx <- dnbinom(0:max(bg_peaks), size=exp_fit$estimate[1], mu=exp_fit$estimate[2])
    if (gene_plot == T) {
      hist(bg_peaks, main='Background peak distribution and negative binomial fit', xlab='# of peaks')
      par(new=T)
      plot(0:max(bg_peaks), hx, ylab = "", xaxt='n', yaxt='n', xlab="", col='firebrick4')
    }
    #Probability that X is larger than number of peaks seen here
    gene_pval = pnbinom(length(gene_inds), size=exp_fit$estimate[1], mu=exp_fit$estimate[2], lower.tail=F)
    gene_peaks[i] = length(gene_inds)
    gene_pvals[i] = gene_pval
    totals[i] = nrow(peak_data)
    i = i + 1
  }
  return(list(gene_peaks, gene_pvals))
}

#Grab all current peak files
peak_files = dir(path=paste(getwd(), '/dnase_data/', sep=''), pattern='*.chr2.sorted.bed.nochr.bed')

#Includes the gene and upstream genomic region reported in Shimizu, et al., so no 1000kb region added
abca_start = 214931541
abca_end = 216017439
#To include potential regulatory regions, add 1000kb
znf_start = 184598366 - 1000000
znf_end = 184939492 + 1000000

abca_pe = get_peak_enrichment(peak_files, abca_start, abca_end, FALSE)
znf_pe = get_peak_enrichment(peak_files, znf_start, znf_end, FALSE)

###Overlap of peaks and variants analysis

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
  #successes
  x = overlap_peaks
  #total successes in background (assume half probability for variant to intersect peak)
  q = num_vars / 2
  #total failures in background
  p = peak_total - overlap_peaks
  #sample size
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
overlap_files = dir(path=paste(getwd(), '/dnase_data/', sep=''), pattern='*.intersected.bed')
abca_overlaps = rep(0, length(overlap_files))
total_overlaps = rep(0, length(overlap_files))
i = 1
for (inFile in overlap_files) {
  overlap_data = fread(paste('dnase_data/',inFile, sep=''))
  starts = overlap_data[,2]
  ends = overlap_data[,3]
  abca_inds = which(overlap_data[,2] >= abca_start & overlap_data[,3] <= abca_end)
  abca_overlaps[i] = length(abca_inds)
  totals[i] = nrow(overlap_data)
  i = i + 1
}


