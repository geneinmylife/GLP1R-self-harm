#####Steiger filtering#####
library(TwoSampleMR)
exposuredata<-read.csv("/Users/tq20202/Desktop/ins.csv")
outlist<-read.csv("/Users/tq20202/Desktop/outlist.csv")
tmp <- data.frame(
  SNP = exposuredata$SNP,
  beta = exposuredata$beta,
  se = exposuredata$se,
  effect_allele = exposuredata$effect_allele,
  other_allele = exposuredata$other_allele,
  eaf = exposuredata$effect_allele_freq,
  phenotype_id = paste(exposuredata$Gene.protein,exposuredata$Type,sep="_"),
  pval = exposuredata$pval,
  samplesize = exposuredata$samplesize
)
exposure_dat <- format_data(tmp, type="exposure", phenotype_col = "phenotype_id",eaf_col = "eaf",samplesize_col = "samplesize")
exposure_dat<-exposure_dat[exposure_dat$exposure=="GLP1R_pQTL",]
outcomeid=outlist[outlist$Outcome.name=="Hospital treatment following self-harm","ID.study"]
outcome_samplesize=as.numeric(outlist[outlist$Outcome.name=="Hospital treatment following self-harm","N_total"])
ncase=as.numeric(outlist[outlist$Outcome.name=="Hospital treatment following self-harm","N_cases"])
ncontrol=as.numeric(outlist[outlist$Outcome.name=="Hospital treatment following self-harm","N_controls"])
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes= outcomeid)
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
###steiger filtering
testres <- generate_odds_ratios(testres)
testres <- cbind(SNP=1,testres,beta.exp=1,beta.out=1,se.exp=1,se.out=1,eaf.exp=1,eaf.out=1,p.exp=1,p.out=1,samplesize.exposure=1)
for (r in 1: nrow(testres)){
  testres[r,1] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],"SNP"]
  testres[r,(ncol(testres)-8):ncol(testres)] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],c("beta.exposure","beta.outcome","se.exposure",
                                                                                                          "se.outcome", "eaf.exposure","eaf.outcome","pval.exposure","pval.outcome",
                                                                                                          "samplesize.exposure")]
}
#testres <- na.omit(testres)
testres <- cbind(testres,rsq.exposure=1,rsq.outcome=1,samplesize.outcome=outcome_samplesize)
testres$rsq.exposure <- (get_r_from_pn(p=testres$p.exp, n=testres$samplesize.exposure))^2
testres$rsq.outcome <- (get_r_from_lor(lor=log(testres$or),af=testres$eaf.out,ncase=ncase,ncontrol=ncontrol,prevalence=ncase/outcome_samplesize))^2
st <- psych::r.test(
  n = testres$samplesize.exposure, 
  n2 = outcome_samplesize, 
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
harmname<-paste("harmdata",outcomeid,sep="_")
harmname<-paste(harmname,".txt",sep="")
write.table(harmdat, file=harmname, row.names=F, col.names=T, sep="\t", quote=F)




