#Toy pipeline to understand how to visualize differences between multiple sequences

require(Biostrings)
require(msa)
library(graph)
library(markovchain)
library(shape)

#Calculate the transition matrix for ACTG- based on all sequences
calcTransitionMat <- function(aligns) {
  nuc_count = matrix(rep(0,5), nrow=1, ncol=5)
  colnames(nuc_count) = c('A', 'T', 'C', 'G', '-')
  numSeqs = length(allSeqs)
  countMat = matrix(rep(0,25), nrow=5, ncol=5)
  colnames(countMat) = c('A', 'T', 'C', 'G', '-')
  rownames(countMat) = c('A', 'T', 'C', 'G', '-')
  for (i in 1:length(aligns)) {
    align = unlist(aligns[i])
    for (j in 1:(length(align)-1)) {
      current_nuc = as.character(align[j])
      next_nuc = as.character(align[j+1])
      nuc_count[,current_nuc] = nuc_count[,current_nuc] + 1
      countMat[current_nuc, next_nuc] = countMat[current_nuc, next_nuc] + 1
    }
  }
  countMat = countMat / as.numeric(nuc_count)
  return(countMat)
}

#Generate toy sequences of length 10
sn = 10
seq1 = sample(DNA_ALPHABET[1:4], size=sn, replace=TRUE)
seq1 = paste(seq1, collapse="")

seq2 = sample(DNA_ALPHABET[1:4], size=sn, replace=TRUE)
seq2 = paste(seq2, collapse="")

seq3 = sample(DNA_ALPHABET[1:4], size=sn, replace=TRUE)
seq3 = paste(seq3, collapse="")

allSeqs = c(seq1, seq2, seq3)
#Align multiple sequences using clustalW
alignment <- msa(allSeqs, "ClustalW", type="dna")

#Calculate transition matrix based on the above sequences
testMat = calcTransitionMat(alignment@unmasked)

al1 = unlist(strsplit(as.character(alignment@unmasked[[1]]), ""))
al2 = unlist(strsplit(as.character(alignment@unmasked[[2]]), ""))
al3 = unlist(strsplit(as.character(alignment@unmasked[[3]]), ""))

#Get consensus alignment
consensus = consensusString(consensusMatrix(alignment))
con_list = unlist(strsplit(consensus, ""))
#Get path and weights based on consensus sequence
#Whenever consensus differs (i.e. = '?'), path splits, and weight splits accordingly
weights = c()
path = c()
for (i in 1:length(con_list)) {
  con_nuc = con_list[i]
  if (con_nuc != '?') {
    path = c(path, con_nuc)
    weights = c(weights, 1)
  }
  else {
    all_nucs = c(al1[i], al2[i], al3[i])
    path = c(path, list(unique(all_nucs)))
    counts = table(all_nucs)
    new_weights = as.data.frame(counts / sum(counts))$Freq
    weights = c(weights, list(new_weights))
  }
}

#Collapse single nucleotide nodes for efficiency
#If several nodes in a row have no divergence, combine
#i.e. 'A' -> 'G' -> 'T' becomes one node of 'AGT'
i = 1
collapsed_nodes = c()
new_weights = c()
while (T) {
  current_node = path[[i]]
  if (length(current_node) > 1) { 
    collapsed_nodes = c(collapsed_nodes, list(path[[i]]))
    new_weights = c(new_weights, weights[i])
    if (i < length(path)) { i = i + 1 }
    else { break }
  }
  else { 
    next_node = c(current_node) 
    single_node = NULL
    while ((length(current_node) == 1) & (i < length(path))) {
      i = i + 1
      current_node = path[[i]]
      if (length(current_node) == 1) { next_node = c(next_node, current_node) }
      else { single_node = path[[i]] }
    }
    collapsed_nodes = c(collapsed_nodes, paste(next_node, collapse=""))
    new_weights = c(new_weights, 1)
    if (!is.null(single_node)) { 
      collapsed_nodes = c(collapsed_nodes, list(single_node)) 
      new_weights = c(new_weights, weights[i])
    }
    if (i >= length(path)) { break }
    else { i = i + 1 }
  }
}

#Create rough visualization of multiple alignment
#Numbers hardcoded currently, may not adapt well to longer/different sequences
emptyplot(c(0, 150))
j = 2
for (i in 1:length(collapsed_nodes)) {
  node = collapsed_nodes[[i]]
  if (length(node) == 1) {
    plotcircle(mid = c(j, 25), r = 20, col = "white")
    textflag(mid = c(j, 25), radx = 0.75, rady = 0.1, lab = node, 
             cex = 0.6, lcol="white", col=NULL)
    j = j + 45
  }
  else {
    h = 25
    for (k in 1:length(node)) {
      plotcircle(mid = c(j, h), r = 20, col = "white")
      textflag(mid = c(j, h), radx = 0.75, rady = 0.1, lab = node[k], 
               cex = 0.6, lcol="white", col=NULL)
      h = h + 45
    }
    j = j + 45
  }
}



