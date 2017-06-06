### Margaret Antonio
### GENE245 Project with Nicole Ferraro

### Functions to (a) Read genotype file as data frame for two classes
###              (b) Train KNN model 
###              (c) Output ROC curves
### Used for Fig.1 Panels A and B


# Packages
library(caret)
library(e1071)
library(pROC)

# Function to format a tab delim genotype file (col=samples, row=rsids)
# Call on input file. Only keeps individuals of pop1 and pop2
# INPUT genoFile: use get_1000genomes_data.sh to convert vcf to genotype tab delim file
# INPUT pop1,pop2: 1000 genomes population codes, {"EAS","SEA","AFR","EUR","AMR"}
# OUTPUT data frame with rows as samples, columns as rsids, cells as numeric genotypes

getPops <- function(genoFile, pop1, pop2){

  tgp.genos <- read.table(genoFile, header = TRUE)
  info<-read.csv("~/Documents/stanford/classes/spring17/bmi245/project/graph_genome_GENE245/data/1000genomes_samples_clean.csv")
  key <- data.frame(sample = info$Sample.name, superpop = info$Superpopulation.code)
  
  # Invert matrix so individuals are rows and columns are variable features (rsids)
  tgp.genos.mat <- as.matrix(tgp.genos[,2:ncol(tgp.genos)])
  tgp.num <- apply(tgp.genos.mat,1,as.numeric)
  tgp.df <-as.data.frame(tgp.num)
  colnames(tgp.df)<-tgp.genos[,1]
  rownames(tgp.df)<-colnames(tgp.genos)[2:ncol(tgp.genos)]
  
  # Remove rsids with NAs
  tgp.df.comp <- tgp.df[ , colSums(is.na(tgp.df)) == 0]
  
  # Add population designation to each individual (row)
  tgp.full <- merge(tgp.df.comp, key, by.x=0, by.y="sample")
  rownames(tgp.full) <-tgp.full[,1]
  tgp.full[,1]<-NULL
  
  # Variation for EUR and AFR only
  tgp.full.eur <- tgp.full[which(tgp.full$superpop == pop1),]
  tgp.full.afr <- tgp.full[which(tgp.full$superpop == pop2)  ,]
  tgp.full.2 <- rbind(tgp.full.eur,tgp.full.afr)
  
  # Get rid of factors without records
  tgp.full.2$superpop <- as.character(tgp.full.2$superpop)
  tgp.full.2$superpop <- as.factor(tgp.full.2$superpop)
  
  # Remove variant sites where ALL individuals are the same
  # uniquelength <- sapply(tgp.full.2,function(x) length(unique(x)))
  # tgp.full.2 <- subset(tgp.full.2, select=uniquelength>1)

  return(tgp.full.2)
  
}

# Function to fit KNN model based on genotype data
# INPUT: whole or subset of data frame outputted by getPops()
# OUTPUT: returns list of the KNN model, training set, and test sets
fitKNN <- function(dataMatrix){
  set.seed(400)
  test.ind <- createDataPartition(y = dataMatrix$superpop, p = 0.75, list = FALSE)
  train <- dataMatrix[-test.ind,]
  test <- dataMatrix[test.ind,]
  ctrl <- trainControl(method='repeatedcv',repeats = 3, number = 10, 
                       p = 0.75, classProbs = TRUE, summaryFunction = twoClassSummary, 
                       savePredictions = TRUE)
  knnFit<- train(superpop ~., data = train, method = "knn", 
                 trControl = ctrl, metric = "ROC", tuneLength = 20)
  return(list(knnFit, train, test))
}

# Function: gets ROC for model trained in fitKNN() and tested here
# Also plots ROC
# INPUT: knnFit, model trained in fitKNN()
#         test, test set returned by fitKNN() 
#         title , string, that will be main title for plot
#         color, string, color for ROC line
#         num, numeric, number of line for plot. If >1 then adds line to plot
# OUTPUT: plots ROC and returns ROC info

getROC <- function(knnFit, test, title, lty = 1, color, num){
  knn.probs <- predict(knnFit,newdata = test, type="prob" )
  knn.ROC <- roc(predictor=knn.probs$EAS,
                 response=test$superpop,
                 levels=rev(levels(test$superpop)))
  if (num==1){
    plot(knn.ROC, main = title, col = color, lty = lty, add = FALSE)
  }
  else{
    plot(knn.ROC, col = color, lty = lty, add = TRUE)
  }
  return(knn.ROC)
}

# Function: creates ROC plot with random subsets of variants
# INPUT data is a matrix created by getPops()
# INPUT max is the largest variant set
# INPUT title is the name of the plot
# OUTPUT prints plot, can return ROC curve info

randSubsetROC <- function(data, max = ncol(data)-1, title){
  for (i in 1:10){
    num.snps = c(5, 10, 20, 30, 40, 50, 75, 100, 500, 1000, max)
    colors = c("red","blue","purple","chocolate4","aquamarine","darkmagenta",
               "deepskyblue3","chartreuse3","coral3","darkslategray", "gold4")
    snps.data <- runif(num.snps[i], min=1, max=ncol(data)-1)
    knn.data <- fitKNN(data[,c(snps.data,ncol(data))])
    knnFit.data <- knn.data[[1]]
    train.data <- knn.data[[2]]
    test.data <- knn.data[[3]]
    print(knnFit.data)
    knn.ROC<- getROC(knnFit.data, test = test.data, color = colors[i], num = i, title)
  }
  
  legend("bottomright",fill=c("red","blue","purple","chocolate4","aquamarine","darkmagenta",
                              "deepskyblue3","chartreuse3","coral3","darkslategray", "gold4"), 
         legend=c("5", "10", "20", "30", "40", "50", "75", "100", "500", "1000", "All 5831"),
         title="# of SNPs", box.col="transparent")
}


# Function: adds a special dotted line to the ROC curve plot created in
# randSubsetROC(), which should be called first
# INPUT data is a data frame created by getPops()
# INPUT subset is a numeric vector of the indices of features to include
addVariantSet <- function(data, subset){
  knn.data <- fitKNN(data[,c(subset,ncol(data))])
  knnFit.data <- knn.data[[1]]
  train.data <- knn.data[[2]]
  test.data <- knn.data[[3]]
  print(knnFit.data)
  knn.ROC<- getROC(knnFit.data, test.data, color = "black", 
                   lty = 2, num = 10, title="")
  
}


###### Calling these functions to produce plots

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

# Highly differentiated SNPs
# Based on Colonna et. al. 2014
rs.abca12 = which(colnames(tgp.abca12)=="rs10180970")
rs.znf = c(which(colnames(tgp.znf)=="rs1344706"),which(colnames(tgp.znf)=="rs4667001"))

# What to call the plot
title = "KNN Classification of Africans and Europeans\n by Biallelic Variant Sites"
title.rand.afr.eur = "KNN classification of EUR and AFR populations\n by random SNPs in chromosome 2"
title.abca12.afr.eur = "KNN classification of EUR and AFR populations\n by random SNPs in ABCA12"

##### Figure 1 Panels A and B
# Panel A
randSubsetROC(tgp.rand, ncol(tgp.abca12)-1, title = title.rand.afr.eur)
addVariantSet(tgp.rand,rs.abca12)
# Panel B
randSubsetROC(tgp.abca12, ncol(tgp.abca12)-1, title = title.abca12.afr.eur)
addVariantSet(tgp.abca12,rs.abca12)

