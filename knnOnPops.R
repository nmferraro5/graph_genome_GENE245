
# Load data. RSIDs in row. Individuals in each column. Fields have 0,1,or 2 for genotype
tgp.genos <- read.table("data/ABCA12/ALL.chr2.214931542-215138591.genos.rsids.txt", header = TRUE)

# Invert matrix so individuals are rows and columns are variable features (rsids)
ncol(tgp.genos)
# 5914 rsids
nrow(tgp.genos)
# 2505 individuals
tgp.genos.inv <- t(tgp.genos[,2:ncol(tgp.genos)])
colnames(tgp.genos.inv)<-tgp.genos[,1]
# Make numeric


# Add outcomes (population) in final column (5915 th column has outcomes)

info<-read.csv("data/1000genomes_samples_clean.csv")
colnames(info)[1]<-"sample.id"
key <- data.frame(sample = info$sample.id, superpop = info$Superpopulation.code)
tgp.genos.full <- merge(tgp.genos.inv, key, by.x=0, by.y="sample")

# Variation for EUR and AFR only
tgp.genos.full.eur <- tgp.genos.full[which(tgp.genos.full$superpop == "EUR")  ,]
tgp.genos.full.afr <- tgp.genos.full[which(tgp.genos.full$superpop == "AFR")  ,]
tgp.genos.full <- rbind(tgp.genos.full.eur,tgp.genos.full.afr)




#rownames(tgp.genos.full)<-tgp.genos.full[,1]
#tgp.genos.full[,1]<-NULL

# Take out outcomes for var set
tgp.vars <- as.matrix(tgp.genos.full[,1:ncol(tgp.genos.full)-1])
# Store rsids
rsids<-colnames(tgp.vars)
colnames(tgp.vars)<-NULL
# Make numeric
tgp.vars.matrix <- as.matrix(tgp.vars)
tgp.vars.matrix.num<-apply(tgp.vars.matrix,1,as.numeric)
# Gets inverted so revert
tgp.vars.matrix.num.inv<-t(tgp.vars.matrix.num)
is.numeric(tgp.vars.matrix.num.inv)

# Get rid of NAs in columns. If get rid in rows then goes down to 418 individuals.
tgp.vars.full<-tgp.vars.matrix.num.inv[ , colSums(is.na(tgp.vars.matrix.num.inv)) == 0]

#Create training and test sets

set.seed(123)
test <- sample.int(nrow(tgp.vars.full), 100)
train.tgp <- tgp.vars.full[-test,]
test.tgp <- tgp.vars.full[test,]

# Outcomes (population membership)

train.pop <- factor(tgp.genos.full$superpop[-test])
test.pop <- factor(tgp.genos.full$superpop[test])

# Run KNN

knn.20 <- knn(train.tgp, test.tgp, train.pop, k=20)
table(knn.20 ,test.pop)


# Redo with Cross Validation
library(caret)
set.seed(1)
outcomes = as.character(tgp.genos.full$superpop)
vars = tgp.vars.full[,1:ncol(tgp.vars.full)]
index <- createFolds(outcomes, k=10)

# Check distribution of individuals into folds
# sapply(index, length)

# Labels
# outcomes[idx[[1]]]
# Vars
# vars[idx[[1]],]
sapply(index, function(i) table(outcomes[i]))

library(rafalib)

for (i in 1:10) {
  pred <- knn(train=vars[ -index[[i]] , ], test=vars[ index[[i]], ], cl=outcomes[ -index[[i]] ], k=40)
  print(paste0(i,") error rate: ", mean(outcomes[ index[[i]] ] != pred)))
}

# Get best k
set.seed(1)
ks <- seq(1,100,2)
res <- sapply(ks, function(k) {
  res.k <- sapply(seq_along(idx), function(i) {
    pred <- knn(train=Xsmall[ -idx[[i]], ],
    test=Xsmall[ index[[i]], ], outcomes[ -index[[i]] ], k = k)
    mean(outcomes[ index[[i]] ] != pred)
    })
  mean(res.k)
})

# Plot k and misclassification
plot(ks, res, pch=16, type ="l",xlab="k in KNN", ylab="Misclassification Error", 
     main = "k-fold Cross Validation, KNN\nAFR and EUR populations from 
     1000 Genomes Project in ABCA12 gene")

# ROC curve (FP rate vs TP rate)

#WORKING ON THIS RIGHT NOW
knn.40 <- knn(train=vars[ -index[[1]] , ], test=vars[ index[[1]], ], outcomes[ -index[[1]] ], k=40,prob=TRUE)
prob <- attr(knn.40, "prob")
