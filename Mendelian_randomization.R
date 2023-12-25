#######GLP1R self-harm code#############
########################################

######MR#######
rm(list=ls())
gc()
library(TwoSampleMR)
library(coloc)
library(vroom)
library(MendelianRandomization)

route<-"/Users/tq20202/Desktop/chrispro/GLP1R_project"
setwd(route)
data<-read.csv("instruments.csv")
data <- data.frame(SNP = data$SNP,
                   chr = data$chr,
                   pos = data$pos38,
                   beta = as.numeric(data$beta),
                   se = as.numeric(data$se),
                   effect_allele = data$effect_allele,
                   other_allele = data$other_allele,
                   eaf=data$eaf,
                   pval = as.numeric(data$pval),
                   Phenotype = data$Phenotype,
                   samplesize = data$samplesize,
                   info = paste(data$Study,data$cis_trans,sep="_"))
data <- format_data(data, type="exposure", phenotype_col = "Phenotype",chr_col = "chr",pos_col = "pos",samplesize_col = "samplesize",info_col = "info")
# LD clumping 
data <-clump_data(data, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1,pop = "EUR")
data <- cbind(data,fstatistics=1)
for (s in 1:nrow(data)){
  z <- data[s,"beta.exposure"]/data[s,"se.exposure"]
  pve <- z^2/(data[s,"samplesize.exposure"]+z^2)
  data[s,"fstatistics"] <- (data[s,"samplesize.exposure"]-2)*pve/(1-pve)
}
print(min(data$fstatistics))
print(max(data$fstatistics))
data<-data[data$fstatistics>10,]

af<-vroom("/Users/tq20202/Desktop/chrispro/brainpro/OneDrive_1_3-28-2023/AF.txt")
f<-list.files("/Users/tq20202/Desktop/chrispro/brainpro/OneDrive_1_3-28-2023/QSM_IDPs_disco")
for (i in 1:length(f)){
  fname<-paste("/Users/tq20202/Desktop/chrispro/brainpro/OneDrive_1_3-28-2023/QSM_IDPs_disco/",f[i],sep="")
  temp<-vroom(fname)
  temp$eaf<-af$af
  outcome<-temp[temp$rsid %in% data$SNP,]
  outcome_samplesize<-19720
  outcome$p<-10^(-1*outcome$`pval(-log10)`)
  outcome<-outcome[!duplicated(outcome$rsid),]
  inter<-intersect(data$SNP,outcome$rsid)
  exposure_dat <- data[data$SNP %in% inter,]
  outcome_dat <- outcome[outcome$rsid %in% inter,]
  outcome_dat <- data.frame(SNP = outcome_dat$rsid,
                            chr = outcome_dat$chr,
                            pos = outcome_dat$pos,
                            beta = as.numeric(outcome_dat$beta),
                            se = as.numeric(outcome_dat$se),
                            effect_allele = outcome_dat$a1,
                            other_allele = outcome_dat$a2,
                            eaf = outcome_dat$eaf,
                            pval = as.numeric(outcome_dat$p),
                            samplesize = outcome_samplesize)
  outcome_dat <- format_data(outcome_dat, type="outcome", chr_col = "chr",pos_col = "pos",samplesize_col = "samplesize")
  harmdat <- harmonise_data(exposure_dat, outcome_dat, action=2)
  test<-harmdat
  test$new<- paste(test$SNP,test$exposure,sep="_")
  test<-test[,-which(colnames(test)=="id.exposure")]
  test <- format_data(test, type="exposure", snp_col = "SNP", pval_col = "pval.exposure",
                      beta_col = "beta.exposure", se_col = "se.exposure",
                      effect_allele_col = "effect_allele.exposure",
                      other_allele_col = "other_allele.exposure",
                      phenotype_col = "new",samplesize_col = "samplesize.exposure")
  testout  <- outcome_dat
  inter<-intersect(test$SNP,testout$SNP)
  test <- test[test$SNP %in% inter,]
  testout <- testout[testout$SNP %in% inter,]
  testharm <- harmonise_data(test, testout, action=2)
  testharm$mr_keep<-TRUE
  testres<-mr(testharm)
  ###steiger filtering
  testres <- generate_odds_ratios(testres)
  testres <- cbind(SNP=1,testres,beta.exp=1,beta.out=1,se.exp=1,se.out=1,p.exp=1,p.out=1,samplesize.exposure=1)
  for (r in 1: nrow(testres)){
    testres[r,1] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],"SNP"]
    testres[r,(ncol(testres)-6):ncol(testres)] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],c("beta.exposure","beta.outcome","se.exposure",
                                                                                                            "se.outcome", "pval.exposure","pval.outcome",
                                                                                                            "samplesize.exposure")]
  }
  #testres <- na.omit(testres)
  testres <- cbind(testres,rsq.exposure=1,rsq.outcome=1,samplesize.outcome=outcome_samplesize)
  testres$rsq.exposure <- (get_r_from_pn(p=testres$p.exp, n=testres$samplesize.exposure))^2
  testres$rsq.outcome <- (get_r_from_pn(p=testres$p.out, n=testres$samplesize.outcome))^2
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
  harmdat<-harmdat[harmdat$steiger_dir==1,]
  type=strsplit(f[i],"[.]")[[1]][1]
  harmname<-paste("harmdata",type,sep="_")
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
  resname<-paste("mrres",type,sep="_")
  resname<-paste(resname,".txt",sep="")
  write.table(mrres, file=resname, row.names=F, col.names=T, sep="\t", quote=F)
}

