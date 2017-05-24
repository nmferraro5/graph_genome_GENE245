install.packages("googleVis")
library(googleVis)
options(gvis.plot.tag='chart')
setwd("~/Documents/stanford/classes/spring17/bmi245/project/graph_genome_GENE245/sankey/")

eur<-read.table("../data/ABCA12/ALL.chr2.214931542-215138591.short.EUR.vcf")
colnames(eur) <- c("chr","pos","id","ref","alt","maf")
eur$source<-"source"
eur$target<-"target"

# Change odd rows to ref 
eur[c(TRUE,FALSE), "source"] <- paste0(eur[c(TRUE,FALSE), "id"],"-",eur[c(TRUE,FALSE), "ref"])
# Change even rows to alt
eur[!c(TRUE,FALSE), "source"] <- paste0(eur[!c(TRUE,FALSE), "id"],"-",eur[!c(TRUE,FALSE), "alt"])

# Should've doen id
eur$id<-eur$source

#Explicitly state maf and ref
eur[c(TRUE,FALSE), "maf"] <- 1-eur[!c(TRUE,FALSE), "maf"]

# Cut out unneeded cols
simple<-eur[,c("source","target","maf")]

# Add target consensus sequence (nodes 1:(nrow(data)/2))
simple$target<-rep(1:(nrow(simple)/2), each = 2)

#Looks like
#> head(simple)
#source target    maf
#1  rs80044567-C      1 0.9871
#2  rs80044567-T      1 0.0129
#3  rs72937917-C      2 0.9940
#4  rs72937917-T      2 0.0060
#5 rs539446331-A      3 1.0000
#6 rs539446331-G      3 0.0000

# Add originating consensus sequence

#Actually rsid->cns should all be freq=1
cns=simple
cns$maf<-1
#> head(cns)
#source target maf
#1  rs80044567-C      1   1
#2  rs80044567-T      1   1
#3  rs72937917-C      2   1
#4  rs72937917-T      2   1

#Fix cns->rs transition matrix

# Swap source and target
simple$temp <-simple$source
simple$source <-simple$target
simple$target <-simple$temp
simple$temp<-NULL

#Decrement source node to capture the cns node before the rs
simple$source<-simple$source - 1

# Add "cns" to each cns node

simple$source<-paste0("cns-",simple$source)
#> head(simple)
#source        target    maf
#1  cns-0  rs80044567-C 0.9871
#2  cns-0  rs80044567-T 0.0129

cns$target <- paste0("cns-",cns$target)
#> head(cns)
#source target maf
#1  rs80044567-C  cns-1   1
#2  rs80044567-T  cns-1   1


# Bind incoming and outgoing rs nodes
pop.eur <- rbind(cns,simple)

