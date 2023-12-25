#######GLP1R self-harm code#############
########################################

###colocization######
rm(list=ls())
gc()
library(TwoSampleMR)
library(coloc)
#coloc function#
coloc.analysis <- function(snp,beta1,beta2,se1,se2,MAF1,MAF2,N1,N2,s){
  dataset1 <- list(snp=snp,beta=beta1, varbeta=se1^2, MAF=MAF1,type="quant", N=N1)
  #type cc (case control) for binary study
  dataset2 <- list(snp=snp,beta=beta2, varbeta=se2^2,MAF=MAF2, type="cc",s=s,N=N2)
  #setting the prior probabilities for association with each trait (p1, p2) and both traits together (p12)
  result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)  
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  names(df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf")
  return(df)
}

res<-read.csv("colocdata.csv")
outcome_id<-c("ieu-a-1126","ukb-b-8777","ieu-a-966","ieu-a-1120","ieu-b-85")
outcome_name<-c("Breast cancer","Cervical cancer","Lung cancer","Ovarian cancer","Prostate cancer")
outcome_samplesize<-c(228951,462933,27209,66450,140254)
ncase<-c(122977,1889,11348,25509,79148)
ncontrol<-c(105974,461044,15861,40941,61106)
bres<-c()
c6<-vroom("/Users/tq20202/reference/GTEX_lookup/look_chr6.csv")
c6$pos37<-sapply(strsplit(c6$variant_id_b37,"_"),"[",2)
for (i in 1: nrow(res)){
  rsid=fdrres[i,"SNP"]$SNP
  chr=fdrres[i,"chr"]$chr
  pos=fdrres[i,"pos37"]$pos37
  outid="ebi-a-GCST006941"
  cc<-c6[c6$pos37>pos-500000 & c6$pos37<pos+500000,]
  exp=extract_outcome_data(snps = cc$rs_id_dbSNP151_GRCh38p7,outcomes= "ubm-a-324")
  out=extract_outcome_data(snps = cc$rs_id_dbSNP151_GRCh38p7,outcomes= "ukb-d-20554_3")
  exp=data.frame(SNP=exp$SNP,effect_allele=exp$effect_allele.outcome,other_allele=exp$other_allele.outcome,effect_allele_freq=exp$eaf.outcome,beta=exp$beta.outcome,se=exp$se.outcome,p=exp$pval.outcome,n=exp$samplesize.outcome)
  out=data.frame(SNP=out$SNP,effect_allele=out$effect_allele.outcome,other_allele=out$other_allele.outcome,effect_allele_freq=out$eaf.outcome,beta=out$beta.outcome,se=out$se.outcome,p=out$pval.outcome,n=out$samplesize.outcome,case=2125)
  temp=intersect(exp$SNP,out$SNP)
  exp<-exp[exp$SNP %in% temp,]
  exp <- exp[order(exp$SNP),]
  exp <- exp[!duplicated(exp$SNP),]
  out<-out[out$SNP %in% temp,]
  out <- out[order(out$SNP),]
  out <- out[!duplicated(out$SNP),]
  df <- data.frame(SNP=exp$SNP,beta1=as.numeric(exp$beta),beta2=as.numeric(out$beta),se1=as.numeric(exp$se),se2=as.numeric(out$se),MAF1=as.numeric(exp$effect_allele_freq),MAF2=as.numeric(out$effect_allele_freq),N1=exp$n,N2=out$n,s=out$case/out$n)
  df <- na.omit(df)
  result <- coloc.analysis(df$SNP,df$beta1, df$beta2, df$se1, df$se2, df$MAF1, df$MAF2, N1=df$N1, N2=df$N2, s=df$s) 
  result <- data.frame(result,num=i,gene=gene,cell=cell,cancer=cancer,cancerid=outid)
  bres <- rbind(bres,result)
}

#pwcoco#
#module load tools/cmake/3.16.3-gcc.9.3
#./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 prepdata/exp.txt  --sum_stats2  prepdata/out.txt  --out res/glp1r_res --out_cond




