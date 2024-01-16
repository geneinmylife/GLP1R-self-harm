#######GLP1R self-harm code#############
########################################

#####ldcheck######
library(TwoSampleMR)
library(MendelianRandomization)
library(ieugwasr)
library(vroom)

setwd("./data/coloc_LDcheck")


rsid="rs880347"
chr="chr6"
pos=39064045
pos37=39031821
f<-vroom("rs880347_500kbregion_pancreas_eqtl_colocdata.csv")
f<-f[f$chr==chr & f$pos38>pos-500000 & f$pos38<pos+500000,]
assoc<-f
  if (nrow(assoc)!=0){
    assoc <- assoc[order(assoc$p),]
    data <- assoc[assoc$p<1E-3,]
    if(nrow(data)>=500){data <- data[1:499,] }
    print(nrow(data))
    
    if (nrow(data)!=0){
      data<-na.omit(data)
      if ("RSID" %in% colnames(data)){
        snp <- append(data$RSID, rsid)
      }else if ("SNP" %in% colnames(data)){
        snp <- append(data$SNP, rsid)
      }
      a <- NULL
      attempts <- 0
      while(attempts<=10){    
        a <- ld_matrix(variants=snp,pop = "EUR")
        if(is.null(a)){attempts<-attempts+1}else{break}
      }
      if(is.null(nrow(a))==TRUE){c<-cbind(snp=rsid, ld_snp=rsid, ld_r2=1)} 
      else {col.index <- which(grepl(rsid,colnames(a)))
      if (length(col.index)>0){
        b <- (a[,col.index])^2
        b <- b[order(b)]
        b <- b[(length(b)-1)] 
        c <- cbind(snp=rsid, ld_snp=names(b), ld_r2=as.numeric(b))
      } else {c <- cbind(snp=rsid, ld_snp="NA", ld_r2="NA")}
      }
    } else {
      c <- c(snp=rsid, ld_snp="NA",ld_r2="NA")}
  } else {c <- c(snp=rsid, ld_snp="NA",ld_r2="NA")}
  





