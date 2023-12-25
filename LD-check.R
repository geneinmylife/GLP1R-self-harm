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




