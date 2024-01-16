#######GLP1R self-harm code#############
########################################

######MR#######
rm(list=ls())
gc()
library(TwoSampleMR)
library(coloc)
library(vroom)
library(MendelianRandomization)

setwd("./data/MR")
data<-read.csv("multi-omic_qtl_instruments.csv")

# construct exposure dataset 
data <- data.frame(SNP = data$SNP,
                   chr = data$chr,
                   pos = data$pos38,
                   beta = as.numeric(data$beta),
                   se = as.numeric(data$se),
                   effect_allele = data$effect_allele,
                   other_allele = data$other_allele,
                   eaf=data$effect_allele_freq,
                   pval = as.numeric(data$pval),
                   Phenotype = data$Phenotype,
                   samplesize = data$samplesize,
                   info = paste(data$Study,data$cis_trans,sep="_"))
data <- format_data(data, type="exposure", phenotype_col = "Phenotype",chr_col = "chr",pos_col = "pos",samplesize_col = "samplesize",info_col = "info")
data <-clump_data(data,clump_kb=10000,clump_r2=0.001,clump_p1=1,clump_p2=1,pop="EUR")

# F statistics>10 to keep relative strong instruments
exposure_dat <- cbind(data,fstatistics=1)
for (s in 1:nrow(exposure_dat)){
  z <- exposure_dat[s,"beta.exposure"]/exposure_dat[s,"se.exposure"]
  pve <- z^2/(exposure_dat[s,"samplesize.exposure"]+z^2)
  exposure_dat[s,"fstatistics"] <- (exposure_dat[s,"samplesize.exposure"]-2)*pve/(1-pve)
}
print(min(exposure_dat$fstatistics))
print(max(exposure_dat$fstatistics))
exposure_dat<-exposure_dat[exposure_dat$fstatistics>10,]

# read already mapped outcome datasets and conduct 2SMR analysis
outcomelist<-read.csv("outcomelist.csv")
for (i in 1:nrow(outcomelist)){
  outcome_dat <- vroom(outcomelist[i,"filename"])
  outcome_dat <- data.frame(SNP = outcome_dat$SNP,
                            chr = outcome_dat$chr,
                            pos = outcome_dat$pos,
                            beta = as.numeric(outcome_dat$beta.outcome),
                            se = as.numeric(outcome_dat$se.outcome),
                            effect_allele = outcome_dat$effect_allele.outcome,
                            other_allele = outcome_dat$other_allele.outcome,
                            eaf = outcome_dat$eaf.outcome,
                            pval = as.numeric(outcome_dat$pval.outcome),
                            samplesize = outcome_dat$samplesize.outcome)
  outcome_dat <- format_data(outcome_dat, type="outcome", chr_col = "chr",pos_col = "pos",samplesize_col = "samplesize")
  harmdat <- harmonise_data(exposure_dat, outcome_dat, action=2)
  test<-harmdat
  test$new<- paste(test$SNP,test$exposure,sep="_")
  test<-test[,-which(colnames(test)=="id.exposure")]
  test <- format_data(test, type="exposure", snp_col = "SNP", pval_col = "pval.exposure",
                      beta_col = "beta.exposure", se_col = "se.exposure",
                      effect_allele_col = "effect_allele.exposure",
                      other_allele_col = "other_allele.exposure",
                      eaf_col = "eaf.exposure",
                      phenotype_col = "new",samplesize_col = "samplesize.exposure")
  testout  <- outcome_dat
  inter<-intersect(test$SNP,testout$SNP)
  test <- test[test$SNP %in% inter,]
  testout <- testout[testout$SNP %in% inter,]
  testharm <- harmonise_data(test, testout, action=2)
  testharm$mr_keep<-TRUE
  testres<-mr(testharm)
  
  ###steiger filtering to check the direction of the association
  testres <- generate_odds_ratios(testres)
  testres <- cbind(SNP=1,testres,beta.exp=1,beta.out=1,se.exp=1,se.out=1,eaf.exp=1,eaf.out=1,p.exp=1,p.out=1,samplesize.exposure=1)
  for (r in 1: nrow(testres)){
    testres[r,1] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],"SNP"]
    testres[r,(ncol(testres)-8):ncol(testres)] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],c("beta.exposure","beta.outcome","se.exposure",
                                                                                                            "se.outcome", "eaf.exposure","eaf.outcome",
                                                                                                            "pval.exposure","pval.outcome",
                                                                                                            "samplesize.exposure")]
  }
  #testres <- na.omit(testres)
  testres <- cbind(testres,rsq.exposure=1,rsq.outcome=1)
  testres$rsq.exposure <- (get_r_from_pn(p=testres$p.exp, n=testres$samplesize.exposure))^2
  if (outcomelist[i,"N_cases"]=="/"){
    samplesize.outcome=as.numeric(outcomelist[i,"N_total"])
    testres$rsq.exposure <- (get_r_from_pn(p=testres$p.out, n=samplesize.outcome))^2
  }else{
    ncase=as.numeric(outcomelist[i,"N_cases"])
    ncontrol=as.numeric(outcomelist[i,"N_controls"])
    samplesize.outcome=as.numeric(outcomelist[i,"N_total"])
    prevalence=ncase/samplesize.outcome
    testres$rsq.outcome <- (get_r_from_lor(lor=log(testres$or),af=testres$eaf.out,ncase=ncase,ncontrol=ncontrol,prevalence=prevalence))^2
  }
  
  st <- psych::r.test(
    n = testres$samplesize.exposure, 
    n2 = testres$samplesize.outcome, 
    r12 = sqrt(testres$rsq.exposure), 
    r34 = sqrt(testres$rsq.outcome)
  )
  testres$steiger_dir <- testres$rsq.exposure > testres$rsq.outcome
  testres$steiger_pval <- pnorm(-abs(st$z)) * 2
  harmdat$match<-paste(harmdat$SNP,harmdat$exposure,sep="_")
  harmdat$steiger_dir <- 0
  harmdat$steiger_pval <- 0
  for (k in 1:nrow(harmdat)){
    harmdat[k,"steiger_dir"] <- testres[testres$exposure==harmdat[k,"match"],"steiger_dir"][1]
    harmdat[k,"steiger_pval"] <- testres[testres$exposure==harmdat[k,"match"],"steiger_pval"][1]
  }
  #harmdat<-harmdat[harmdat$steiger_dir==1,]
  outtype=outcomelist[i,"Outcome.name"]
  harmname<-paste("harmdata",outtype,sep="_")
  harmname<-paste(harmname,".txt",sep="")
  write.table(harmdat, file=harmname, row.names=F, col.names=T, sep="\t", quote=F)
  exp<-unique(harmdat$exposure)
  exp1=exposure_dat
  exp2=NULL
  mrres <- c()
  for (t in 1:length(exp)){
    dat <- harmdat[harmdat$exposure==exp[t],]
    dat<-dat[!is.na(dat$SNP),]
    if (nrow(dat)==1){
      result <- mr_wald_ratio(b_exp= dat$beta.exposure, b_out=dat$beta.outcome, se_exp= dat$se.exposure, se_out= dat$se.outcome)
      result <- data.frame(id=exp[t],method="Wald_ratio",nsnp=result$nsnp,b=result$b,se=result$se,CIlower=NA,CIupper=NA,
                           pval=result$pval,intercept=NA,intercept_se=NA,inter_CIlower=NA,inter_CIupper=NA,
                           intercept_pval=NA,hetero_Q=NA,hetero_pvale=NA)}
    else if (nrow(dat)==2 & all(dat$SNP %in% exp1$SNP)){rho <- ld_matrix(dat$SNP, pop = "EUR")
    result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                      bxse = dat$se.exposure,
                                                      by = dat$beta.outcome,
                                                      byse = dat$se.outcome,
                                                      cor = rho))
    result <- data.frame(id=exp[t],
                         method="IVW",
                         nsnp=result$SNPs,
                         b=result$Estimate,
                         se=result$StdError,
                         CIlower=result$CILower,
                         CIupper=result$CIUpper,
                         pval=result$Pvalue,
                         intercept=NA,
                         intercept_se=NA,
                         inter_CIlower=NA,
                         inter_CIupper=NA,
                         intercept_pval=NA,
                         hetero_Q=result$Heter.Stat[1],
                         hetero_pvale=result$Heter.Stat[2])}
    else if (nrow(dat)==2 & any(dat$SNP %in% exp2$SNP)){result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                                                                          bxse = dat$se.exposure,
                                                                                                          by = dat$beta.outcome,
                                                                                                          byse = dat$se.outcome))
    result <- data.frame(id=exp[t],
                         method="IVW",
                         nsnp=result$SNPs,
                         b=result$Estimate,
                         se=result$StdError,
                         CIlower=result$CILower,
                         CIupper=result$CIUpper,
                         pval=result$Pvalue,
                         intercept=NA,
                         intercept_se=NA,
                         inter_CIlower=NA,
                         inter_CIupper=NA,
                         intercept_pval=NA,
                         hetero_Q=result$Heter.Stat[1],
                         hetero_pvale=result$Heter.Stat[2])}
    else if (nrow(dat)>2 & all(dat$SNP %in% exp1$SNP)){rho <- ld_matrix(dat$SNP, pop = "EUR")
    result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
                                                        bxse = dat$se.exposure,
                                                        by = dat$beta.outcome,
                                                        byse = dat$se.outcome,
                                                        cor = rho))
    result <- data.frame(id=exp[t],
                         method="MR-Egger",
                         nsnp=result$SNPs,
                         b=result$Estimate,
                         se=result$StdError.Est,
                         CIlower=result$CILower.Est,
                         CIupper=result$CIUpper.Est,
                         pval=result$Pvalue.Est,
                         intercept=result$Intercept,
                         intercept_se=result$StdError.Int,
                         inter_CIlower=result$CILower.Int,
                         inter_CIupper=result$CIUpper.Int,
                         intercept_pval=result$Pvalue.Int,
                         hetero_Q=result$Heter.Stat[1],
                         hetero_pvale=result$Heter.Stat[2])
    result2 <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                       bxse = dat$se.exposure,
                                                       by = dat$beta.outcome,
                                                       byse = dat$se.outcome,
                                                       cor = rho))
    result2 <- data.frame(id=exp[t],
                          method="IVW",
                          nsnp=result2$SNPs,
                          b=result2$Estimate,
                          se=result2$StdError,
                          CIlower=result2$CILower,
                          CIupper=result2$CIUpper,
                          pval=result2$Pvalue,
                          intercept=NA,
                          intercept_se=NA,
                          inter_CIlower=NA,
                          inter_CIupper=NA,
                          intercept_pval=NA,
                          hetero_Q=result2$Heter.Stat[1],
                          hetero_pvale=result2$Heter.Stat[2])
    result <- rbind(result,result2)}
    else if (nrow(dat)>2 & any(dat$SNP %in% exp2$SNP)){result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
                                                                                                           bxse = dat$se.exposure,
                                                                                                           by = dat$beta.outcome,
                                                                                                           byse = dat$se.outcome))
    result <- data.frame(id=exp[t],
                         method="MR-Egger",
                         nsnp=result$SNPs,
                         b=result$Estimate,
                         se=result$StdError.Est,
                         CIlower=result$CILower.Est,
                         CIupper=result$CIUpper.Est,
                         pval=result$Pvalue.Est,
                         intercept=result$Intercept,
                         intercept_se=result$StdError.Int,
                         inter_CIlower=result$CILower.Int,
                         inter_CIupper=result$CIUpper.Int,
                         intercept_pval=result$Pvalue.Int,
                         hetero_Q=result$Heter.Stat[1],
                         hetero_pvale=result$Heter.Stat[2])
    result2 <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                       bxse = dat$se.exposure,
                                                       by = dat$beta.outcome,
                                                       byse = dat$se.outcome))
    result2 <- data.frame(id=exp[t],
                          method="IVW",
                          nsnp=result2$SNPs,
                          b=result2$Estimate,
                          se=result2$StdError,
                          CIlower=result2$CILower,
                          CIupper=result2$CIUpper,
                          pval=result2$Pvalue,
                          intercept=NA,
                          intercept_se=NA,
                          inter_CIlower=NA,
                          inter_CIupper=NA,
                          intercept_pval=NA,
                          hetero_Q=result2$Heter.Stat[1],
                          hetero_pvale=result2$Heter.Stat[2])
    result <- rbind(result,result2)}
    mrres <- rbind(mrres,result)
  }
  #####FDR Adjustment for Pvalue######
  mrres<-generate_odds_ratios(mrres)
  mrres$fdr <- p.adjust(mrres$pval, method = "fdr", n = length(mrres$pval))
  resname<-paste("mrres",outtype,sep="_")
  resname<-paste(resname,".txt",sep="")
  write.table(mrres, file=resname, row.names=F, col.names=T, sep="\t", quote=F)
}




