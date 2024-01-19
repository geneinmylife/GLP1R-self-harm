library(TwoSampleMR)
library("readxl")
rm(list=ls(all=TRUE)) 

#self-harm traits as exposures
ids<- c("ukb-d-20480","ukb-d-20554_1","ukb-d-20554_3","ukb-d-20554_5","ukb-d-20554_6")

#Due to lack of instruments with conventional P value cut of 5e-8, a P value threshold of 1e-6 was used to select instruments for self-harm behaviours. 
exposure_dat <- extract_instruments(ids,p1 = 1e-06,p2 = 1e-06)
#exposure_dat

###load outcome from local files
    try(outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      #13085_18_GLP1R_GLP1R.txt.gz was downloaded from Decode (https://www.decode.com/summarydata/) and merged with effect allele freq from assocvariants.annotated.txt.gz
      filename ="./data/13085_18_GLP1R_GLP1R.txt.gz.tab",
      sep = "\t",
      snp_col = "snp",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      eaf_col = "effect_allele_freq",
      pval_col = "p",
      samplesize_col = "n"
      ))

  outcome_dat$outcome <- "Protein levels of GLP1R"

###run MR analyses
  dat <- NULL
  try(dat <- harmonise_data(exposure_dat, outcome_dat))

  mr_results <- NULL
  mr_hetero <- NULL
  mr_pleio <- NULL
  mr_single <- NULL
  try(mr_results <- mr(dat))  
  mr_hetero <- mr_heterogeneity(dat)
  mr_pleio <- mr_pleiotropy_test(dat) 
  try(mr_single <- mr_singlesnp(dat))
  try(mr_loo <- mr_leaveoneout(dat))
  
  exposure <- "Self-harm-GLP1R"

  
