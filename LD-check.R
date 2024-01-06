#######GLP1R self-harm code#############
########################################
#####ldcheck######
library(TwoSampleMR)
library(MendelianRandomization)
library(ieugwasr)
library(vroom)

route<-"/Users/tq20202/Desktop/chrispro/chris0925/mediation"
setwd(route)

fdrres<-data
posfile<-vroom("/Users/tq20202/positions.txt")
ldres<-c()
for (i in 1:nrow(fdrres)){
  rsid=fdrres[i,"SNP"]
  chr=fdrres[i,"chr.outcome"]
  pos=fdrres[i,"pos.outcome"]
  id=fdrres[i,"exposure"]
  outtype=fdrres[i,"outtype"]
  if (outtype<10) {
    outtype<-paste("000",outtype,sep="")
  }else if (outtype>10 & outtype<100) {
    outtype<-paste("00",outtype,sep="")
  }else if (outtype>100 & outtype<1000) {
    outtype<-paste("0",outtype,sep="")}
  outfile<-paste("/Users/tq20202/ddd/",outtype,".txt",sep="")
  f<-vroom(outfile)
  f<-cbind(posfile,f)
  f<-f[f$CHROM==chr & f$POS>pos-500000 & f$POS<pos+500000,]
  f$p<-10^(-1*f$PVAL)
  assoc<-f
  if (nrow(assoc)!=0){
    assoc <- assoc[order(assoc$p),]
    data <- assoc[assoc$p<1E-3,]
    if(nrow(data)>=500){data <- data[1:499,] }
    print(nrow(data))
    
    if (nrow(data)!=0){
      data<-na.omit(data[,1:10])
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
  
  ldres <- rbind(ldres,c)
}




