#######GLP1R self-harm code#############
########################################

######colocalization#######
rm(list=ls())
gc()
library(TwoSampleMR)
library(coloc)
library(vroom)
setwd("./data/coloc_LDcheck")

##coloc function##
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

coloc.abf <- function(dataset1, dataset2, MAF=NULL, 
                      p1=1e-4, p2=1e-4, p12=1e-5) {
  
  if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
    dataset1$MAF <- MAF
  if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
    dataset2$MAF <- MAF
  
  df1 <- process.dataset(d=dataset1, suffix="df1")
  df2 <- process.dataset(d=dataset2, suffix="df2")
  p1=adjust_prior(p1,nrow(df1),"1")
  p2=adjust_prior(p2,nrow(df2),"2")
  
  merged.df <- merge(df1,df2)
  p12=adjust_prior(p12,nrow(merged.df),"12")
  
  if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")
  
  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)
  ############################## 
  pp.abf <- combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)  
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)
  
  output<-list(summary=results,
               results=merged.df,
               priors=c(p1=p1,p2=p2,p12=p12))
  class(output) <- c("coloc_abf",class(output))
  return(output)
}

process.dataset <- function(d, suffix) {
  nd <- names(d)
  if("beta" %in% nd && "varbeta" %in% nd) {  
    if(d$type=="quant" && !('sdY' %in% nd)) 
      d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
    df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
                              V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    df$snp <- as.character(d$snp)
    if("position" %in% nd)
      df <- cbind(df,position=d$position)
    return(df)
  }
  
  if("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) { 
    df <- data.frame(pvalues = d$pvalues,
                     MAF = d$MAF,
                     N=d$N,
                     snp=as.character(d$snp))    
    snp.index <- which(colnames(df)=="snp")
    colnames(df)[-snp.index] <- paste(colnames(df)[-snp.index], suffix, sep=".")
    abf <- approx.bf.p(p=df$pvalues, f=df$MAF, type=d$type, N=df$N, s=d$s, suffix=suffix)
    df <- cbind(df, abf)
    if("position" %in% nd)
      df <- cbind(df,position=d$position)
    return(df)  
  }
  
  stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(pvalues, MAF, N, type)")
}

combine.abf <- function(l1, l2, p1, p2, p12, quiet=FALSE) {
  stopifnot(length(l1)==length(l2))
  lsum <- l1 + l2
  lH0.abf <- 0
  lH1.abf <- log(p1) + logsum(l1)
  lH2.abf <- log(p2) + logsum(l2)
  lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
  lH4.abf <- log(p12) + logsum(lsum)
  
  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
  my.denom.log.abf <- logsum(all.abf)
  pp.abf <- exp(all.abf - my.denom.log.abf)
  names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
  if(!quiet) {
    print(signif(pp.abf,3))
    print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  }
  return(pp.abf)
}

adjust_prior=function(p,nsnps,suffix="") {
  if(nsnps * p >= 1) { ## for very large regions
    warning(paste0("p",suffix," * nsnps >= 1, setting p",suffix,"=1/(nsnps + 1)"))
    1/(nsnps + 1)
  } else {
    p
  }
}

logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}
logdiff <- function(x,y) {
  my.max <- max(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  return(my.res)
}
sdY.est <- function(vbeta, maf, n) {
  warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[['oneover']]
  if(cf < 0)
    stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
  return(sqrt(cf))
}

approx.bf.estimates <- function (z, V, type, suffix=NULL, sdY=1) {
  sd.prior <- if (type == "quant") { 0.15*sdY } else { 0.2 }
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}


##coloc analysis##
snp="rs877446"
chr="chr6"
pos=39063263
pos37=39031039
outcomeid="ukb-d-20554_1"
outcome_samplesize=117733
ncase=2125
f<-vroom("rs877446_500kbregion_eqtl_T2D_colocdata.csv")
f<-f[f$pos>pos37-500000 & f$pos<pos37+500000, ]
exp=data.frame(SNP=f$SNP,effect_allele=f$effect_allele,effect_allele_freq=f$effect_allele_freq,other_allele=f$other_allele,beta=f$beta,se=f$se,p=f$p,n=f$n)
out <- vroom("rs877446_500kbregion_eqtl_selfharm_ukb-d-20554_1_colocdata.csv")
temp=intersect(exp$SNP,out$SNP)
exp<-exp[exp$SNP %in% temp,]
exp <- exp[order(exp$SNP),]
exp <- exp[!duplicated(exp$SNP),]
out<-out[out$SNP %in% temp,]
out <- out[order(out$SNP),]
out <- out[!duplicated(out$SNP),]
out=data.frame(SNP=out$SNP,effect_allele=out$effect_allele,other_allele=out$other_allele,effect_allele_freq=out$effect_allele_freq,beta=out$beta,se=out$se,p=out$p,n=outcome_samplesize,case=ncase)
exp=data.frame(SNP=exp$SNP,effect_allele=exp$effect_allele,other_allele=exp$other_allele,effect_allele_freq=exp$effect_allele_freq,beta=exp$beta,se=exp$se,p=exp$p,n=exp$n)
df <- data.frame(SNP=exp$SNP,beta1=as.numeric(exp$beta),beta2=as.numeric(out$beta),se1=as.numeric(exp$se),se2=as.numeric(out$se),MAF1=as.numeric(exp$effect_allele_freq),MAF2=as.numeric(out$effect_allele_freq),N1=exp$n,N2=out$n,s=out$case/out$n)
df <- na.omit(df)
result <- coloc.analysis(df$SNP,df$beta1, df$beta2, df$se1, df$se2, df$MAF1, df$MAF2, N1=df$N1, N2=df$N2, s=df$s) 


f<-vroom("rs877446_500kbregion_eqtl_BMI_ieu-b-40_colocdata.csv")
f<-f[f$pos>pos37-500000 & f$pos<pos37+500000, ]
exp=data.frame(SNP=f$SNP,effect_allele=f$effect_allele,effect_allele_freq=f$effect_allele_freq,other_allele=f$other_allele,beta=f$beta,se=f$se,p=f$p,n=f$n)
out <- vroom("rs877446_500kbregion_eqtl_selfharm_ukb-d-20554_1_colocdata.csv")
temp=intersect(exp$SNP,out$SNP)
exp<-exp[exp$SNP %in% temp,]
exp <- exp[order(exp$SNP),]
exp <- exp[!duplicated(exp$SNP),]
out<-out[out$SNP %in% temp,]
out <- out[order(out$SNP),]
out <- out[!duplicated(out$SNP),]
out=data.frame(SNP=out$SNP,effect_allele=out$effect_allele,other_allele=out$other_allele,effect_allele_freq=out$effect_allele_freq,beta=out$beta,se=out$se,p=out$p,n=outcome_samplesize,case=ncase)
exp=data.frame(SNP=exp$SNP,effect_allele=exp$effect_allele,other_allele=exp$other_allele,effect_allele_freq=exp$effect_allele_freq,beta=exp$beta,se=exp$se,p=exp$p,n=exp$n)
df <- data.frame(SNP=exp$SNP,beta1=as.numeric(exp$beta),beta2=as.numeric(out$beta),se1=as.numeric(exp$se),se2=as.numeric(out$se),MAF1=as.numeric(exp$effect_allele_freq),MAF2=as.numeric(out$effect_allele_freq),N1=exp$n,N2=out$n,s=out$case/out$n)
df <- na.omit(df)
result <- coloc.analysis(df$SNP,df$beta1, df$beta2, df$se1, df$se2, df$MAF1, df$MAF2, N1=df$N1, N2=df$N2, s=df$s) 









