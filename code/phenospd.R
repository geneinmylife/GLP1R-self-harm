#######GLP1R self-harm code#############
########################################

#####phenospd######
library(vroom)

setwd("./data/phenospd")

## Read in correlation matrix:  
phenocorr<-vroom("suicide:selfharm_traits.txt",col_names = F)
# Set NA value as 0 to run the SpD analysis
phenocorr[is.na(phenocorr)] <- 0

# For multiple test correction the sign of the correlation is irrelevant (i.e., so we're best to input absolute values)
corr.matrix<-abs(phenocorr)  

## Remove Duplicate Columns:
corr.matrix.RemoveDupCol <- corr.matrix[!duplicated((corr.matrix))]

## Remove Duplicate Rows:
corr.matrix.RemoveDupRow <- unique((corr.matrix.RemoveDupCol))

## Remove Redundant VAR Names:
VARnames.NonRedundant<-as.matrix(dimnames(corr.matrix.RemoveDupCol)[[2]])
colnames(VARnames.NonRedundant)<-"VAR"

evals<-eigen(t(corr.matrix.RemoveDupRow),symmetric=T)$values

oldV<-var(evals)
M<-length(evals)
L<-(M-1)
Meffold<-M*(1-(L*oldV/M^2))

labelevals<-array(dim=M)
for(col in 1:M) { labelevals[col]<-c(col) }
levals<-cbind(labelevals, evals)

newevals<-evals
for(i in 1:length(newevals)) { 
  if(newevals[i] < 0) { 
    newevals[i] <- 0
  }
}

newlevals<-cbind(labelevals, newevals)

newV<-var(newevals)
Meffnew<-M*(1-(L*newV/M^2))

NewResulttemp<-c('Effective Number of Independent Variables [Veff] (-ve values set to zero):', round(Meffnew,dig=4))
NewResulttemp

## Implement improved approach of Li and Ji. Heredity 2005 95:221-227

IntLinewevals<-newevals

for(i in 1:length(IntLinewevals)) {
  if(IntLinewevals[i] >= 1 ) {
    IntLinewevals[i] <- 1
  }
  if(IntLinewevals[i] < 1 ) {
    IntLinewevals[i] <- 0
  }
}

NonIntLinewevals <- newevals-floor(newevals)

MeffLi <- sum(NonIntLinewevals+IntLinewevals)

NewResultLitemp1<-c('Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005):')
NewResultLitemp2<-round(MeffLi,dig=4)
NewResultLitemp1
NewResultLitemp2









