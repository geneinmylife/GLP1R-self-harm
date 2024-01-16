#######GLP1R self-harm code#############
########################################

#######hyprcoloc code######
rm(list=ls())
gc()
library(gwasglue2)
library(VariantAnnotation)
library(gwasvcf)
library(ieugwasr)
library(susieR) 
library(hyprcoloc)
library(genetics.binaRies)
plink_path <- get_plink_binary()
file.exists(plink_path)

###index snp and region######
snp="rs877446"
chr="chr6"
pos=39063263
pos37=39031039

snp="rs880347"
chr="chr6"
pos=39064045
pos37=39031821

snp="rs10305439"
chr="chr6"
pos=39056940
pos37=39024716


setwd("./data/coloc_LDcheck/moloc")
e<-list.files("eqtl")
s<-list.files("sqtl")
bed_ref <- "./data/coloc_LDcheck/moloc/ref/EUR"


sumseteqtl <- lapply(seq_along(e), function(i){
  # create summarysets
  t1<-read.table(file=paste("./eqtl/",e[i],sep=""), header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
  t1<-t1[t1$chr==chr & t1$position>pos-500000 & t1$position<pos+500000,]
  t1<-dplyr:: as_tibble(t1)
  meta1 <- create_metadata(sample_size = 575, trait = "eqtl",id=strsplit(e[i],split="[.]")[[1]][1])
  dt <- create_summaryset(t1, metadata=meta1,build="GRCh38")
  dt<- liftover(dt,to = "GRCh37")
})


sumsetsqtl <- lapply(seq_along(s), function(i){
  # create summarysets
  t1<-read.table(file=paste("./sqtl/",s[i],sep=""), header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
  t1<-t1[t1$chr==chr & t1$position>pos-500000 & t1$position<pos+500000,]
  t1<-dplyr:: as_tibble(t1)
  meta1 <- create_metadata(sample_size = 305, trait = "sqtl",id=strsplit(s[i],split="[.]")[[1]][1])
  dt <- create_summaryset(t1, metadata=meta1,build="GRCh38")
  dt<- liftover(dt,to = "GRCh37")
})


l<-pos37-500000
u<-pos37+500000
region<-paste("6:",l,"-",u,sep="")
outcome <- ieugwasr::associations(variants = region, id = "ukb-d-20480")
sumsetout <- create_summaryset(outcome,qc = TRUE)

summarysets <- c(sumseteqtl, sumsetsqtl, sumsetout)
dataset <-  create_dataset(summarysets, harmonise = TRUE, tolerance = 0.08, action = 1)%>%
  # harmonise dataset against LD matrix
  harmonise_ld(., bfile = bed_ref , plink_bin = plink_path)

res_dataset<-hyprcoloc(dataset = dataset)
print(res_dataset)



