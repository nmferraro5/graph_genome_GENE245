library(data.table)
library(fitdistrplus)
library(ggplot2)

setwd('/Users/nicoleferraro/Documents/Stanford/Spring_1617/STATS345/Project/graph_genome_GENE245/')

###Peak enrichment analysis

#Function to calculate peak enrichment in a gene region as compared to background
#Inputs: list of sample peak files, tissue, gene start and end location, and a boolean if plots are created
#Returns a list of p-vals describing significance of peak enrichment for that region and peak counts
get_peak_enrichment <- function(peak_files, tissue, gene_start, gene_end, gene_plot)
{
  gene_peaks = rep(0, length(peak_files))
  gene_pvals = rep(0, length(peak_files))
  totals = rep(0, length(peak_files))
  gene_length = gene_end - gene_start
  i = 1
  #Consider parallelizing if number of samples increases drastically in future
  for (inFile in peak_files) {
    input_path = paste(paste('dnase_data/', tissue, sep=''), '/', sep='')
    peak_data = fread(paste(input_path,inFile, sep=''))
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
      plot(0:max(bg_peaks), hx, ylab = "", xaxt='n', yaxt='n', xlab="", col='slategrey')
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

#Includes the gene and upstream genomic region reported in Shimizu, et al., so no 1000kb region added
abca_start = 214931541
abca_end = 216017439
#To include potential regulatory regions, add 1000kb
znf_start = 184598366 - 1000000
znf_end = 184939492 + 1000000

peak_files_skin = dir(path=paste(getwd(), '/dnase_data/Skin/', sep=''), pattern='*.chr2.sorted.bed.nochr.bed')

abca_pe_skin = get_peak_enrichment(peak_files_skin, 'Skin', abca_start, abca_end, TRUE)
abca_pe_pvals = abca_pe_skin[[2]]
znf_pe_skin = get_peak_enrichment(peak_files_skin, 'Skin', znf_start, znf_end, FALSE)
znf_pe_pvals = znf_pe_skin[[2]]

peak_files_brain = dir(path=paste(getwd(), '/dnase_data/Brain/', sep=''), pattern='*.chr2.sorted.bed.nochr.bed')

abca_pe_brain = get_peak_enrichment(peak_files_brain, 'Brain', abca_start, abca_end, FALSE)
abca_pe_pvals = abca_pe_brain[[2]]
znf_pe_brain = get_peak_enrichment(peak_files_brain, 'Brain', znf_start, znf_end, FALSE)
znf_pe_pvals = znf_pe_brain[[2]]

###Overlap of peaks and variants analysis

#Read in sample information file
metadata = fread('dnase_metadata_2016-12-05.txt')

#Peak overlap file with summary of all peak-variant overlaps
overlap_skin_data = fread('dnase_data/Skin/peak_overlap.txt')
colnames(overlap_skin_data) = c('Sample', 'Overlap', 'Total')
overlap_skin_data$Proportion = overlap_skin_data$Overlap / overlap_skin_data$Total
sample_ids = gsub("_", ";", overlap_data$Sample)
meta = metadata[which(metadata$encode_expt_accession %in% sample_ids),]
meta_order = meta[order(meta$encode_expt_accession),]

#Look at proportion of overlaps occurring in gene region of interest
overlap_files_skin = dir(path=paste(getwd(), '/dnase_data/Skin/', sep=""), pattern='*.intersected.bed', full.names=T)
overlap_files_brain = dir(path=paste(getwd(), '/dnase_data/Brain/', sep=""), pattern='*.intersected.bed', full.names=T)
overlap_files = overlap_files_brain

#Store variant counts, in peak regions, and in each gene region
avg_length = rep(0, length(overlap_files_brain))
var_peak_all = c()
abca_overlaps = rep(0, length(overlap_files))
znf_overlaps = rep(0, length(overlap_files))

i = 1
#Get count of variants in peak regions in all samples
for (inFile in overlap_files) {
  overlap_data = fread(inFile)
  peaks = overlap_data[,8]
  avg_length[i] = mean(overlap_data$V7 - overlap_data$V6)
  uniq_peaks = unlist(unique(peaks))
  rand_inds = sample(1:length(uniq_peaks), 100)
  var_peak_overlap = unlist(lapply(rand_inds, function(x) length(which(unlist(peaks) == uniq_peaks[x]))))
  var_peak_all = c(var_peak_all, var_peak_overlap)
  #check abca and znf specifically
  abca_inds = which(overlap_data$V6 >= abca_start & overlap_data$V7 <= abca_end)
  znf_inds = which(overlap_data$V6 >= znf_start & overlap_data$V7 <= znf_end)
  abca_overlaps[i] = length(abca_inds)
  znf_overlaps[i] = length(znf_inds)
  i = i + 1
}

#Variant information
variant_data = fread('data/ALL.chr2.sorted.position.vcf')
colnames(variant_data) = c('Chr', 'Start', 'End', 'RSID')
num_vars = nrow(variant_data)

#Sample number of variants in random regions of length x
mean_peak = mean(avg_length)
rand_var_inds = sample(1:max(variant_data$End) - mean_peak, 1600)
num_vars_in_sample = unlist(lapply(rand_var_inds, function(x) length(which(variant_data$Start >= x & variant_data$End <= x+mean_peak))))

hist(num_vars_in_sample, col='royalblue4', xlim=c(0,100), ylim=c(0,700), main='Number of variants in sampled regions - Brain', xlab='Number variants overlapping sample', cex.lab=1.5, cex.main=1.5)
par(new=TRUE)
hist(var_peak_all, col=adjustcolor("firebrick3", alpha.f = 0.6), xlim=c(0,100), ylim=c(0,700), main="", xlab="", breaks=20, ylab='')
legend(50,500, c('Background', 'Peaks'), lty=c(1,1), lwd=c(2,2), col=c("royalblue4",adjustcolor("firebrick3", alpha.f = 0.6)), bty='n', cex=1.5)

#Compare means of two distributions
#Check for equal variance
ftest = var.test(num_vars_in_sample, var_peak_all)
#T test
ttest = t.test(var_peak_all, num_vars_in_sample, var.equal=F, alternative="less")

#Look at abca and znf overlaps
#First fit negative binomial to background distribution
exp_peak_fit = fitdist(num_vars_in_sample, "nbinom")
#Probability that X is larger than number of peaks seen here
gene_pval = pnbinom(length(gene_inds), size=exp_peak_fit$estimate[1], mu=exp_peak_fit$estimate[2], lower.tail=F)

