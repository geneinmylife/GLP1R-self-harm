#######GLP1R self-harm code#############
########################################

#### PGS######
# Clear the work environment
rm(list = ls())

#input data
ivGLP1R<-read.csv("/path_to_dataset/")
phenotype<-read.csv("/path_to_dataset/")
phenotype_other<-read.csv("/path_to_dataset/")
exclude_all<-read.csv("/path_to_dataset/")
exclude_withdraw<-read.csv("/path_to_dataset/")
t2d<-read.csv("/path_to_dataset/")

information<-read.csv("/path_to_dataset/")

#clean data
colnames(phenotype)[1]<-"app"
colnames(phenotype_other)<-c("app","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","sex","height","BMI","weight","age","array","geneticSEX","SBPm","DBPm","DBPa","SBPa")
phenotype_other<-phenotype_other[,c(1:12,14,16:18)]
phenotype_other$array2[phenotype_other$array>=-11 & phenotype_other$array<0]<-1
phenotype_other$array2[phenotype_other$array>0 & phenotype_other$array<=95]<-2
exclude_all_use<-exclude_all[,c("ieu","app")]
colnames(exclude_withdraw)<-"app"
colnames(t2d)[1]<-"app"
t2d$T2D[is.na(t2d$T2D) & is.na(t2d$T1D)]<-0

#exclude participants
GLP1R_use<-merge(ivGLP1R,phenotype,by="app",all.x=TRUE)
GLP1R_use<-merge(GLP1R_use,phenotype_other,by="app",all.x=TRUE)
GLP1R_use<-merge(GLP1R_use,t2d,by="app",all.x=TRUE)
GLP1R_use$check[GLP1R_use$app%in%exclude_all_use$app]<-1
GLP1R_use$check2[GLP1R_use$app%in%exclude_withdraw$app]<-1
GLP1R_use2<-GLP1R_use[which(is.na(GLP1R_use$check) & is.na(GLP1R_use$check2)),]

#check whether non-EUR participants were excluded
exclude_EUR<-read.csv("/path_to_dataset/")
GLP1R_use2_check<-merge(GLP1R_use2,exclude_EUR,by="ieu")
################################################################################
#outcome definition
#ever self-harm
GLP1R_use2$self_harm[GLP1R_use2$X20480.0.0==1]<-1
GLP1R_use2$self_harm[GLP1R_use2$X20480.0.0==0]<-0
table(GLP1R_use2$self_harm,exclude=NULL)

#suicide
GLP1R_use2$suicide[GLP1R_use2$X20480.0.0==0]<-0
GLP1R_use2$suicide[GLP1R_use2$X20483.0.0==1]<-1
table(GLP1R_use2$suicide,exclude=NULL)

#mental health sevice
GLP1R_use2$MHS[GLP1R_use2$X20480.0.0==0]<-0
GLP1R_use2$MHS[GLP1R_use2$X20554.0.1==1]<-1
table(GLP1R_use2$MHS,exclude=NULL)

#hospital
GLP1R_use2$hospital[GLP1R_use2$X20480.0.0==0]<-0
GLP1R_use2$hospital[GLP1R_use2$X20554.0.1==3 | GLP1R_use2$X20554.0.2==3]<-1
table(GLP1R_use2$hospital,exclude=NULL)

#GP
GLP1R_use2$GP[GLP1R_use2$X20480.0.0==0]<-0
GLP1R_use2$GP[GLP1R_use2$X20554.0.1==5 | GLP1R_use2$X20554.0.2==5 | GLP1R_use2$X20554.0.3==5 | GLP1R_use2$X20554.0.4==5]<-1
table(GLP1R_use2$GP,exclude=NULL)

#help
GLP1R_use2$help[GLP1R_use2$X20480.0.0==0]<-0
GLP1R_use2$help[GLP1R_use2$X20554.0.1==6| GLP1R_use2$X20554.0.2==6 | GLP1R_use2$X20554.0.3==6 | GLP1R_use2$X20554.0.4==6 | GLP1R_use2$X20554.0.5==6]<-1
table(GLP1R_use2$help,exclude=NULL)

#BMI
a<-mean(GLP1R_use2$BMI,na.rm=TRUE)
b<-sd(GLP1R_use2$BMI,na.rm=TRUE)
bmi1<-a-4*b
bmi2<-a+4*b
GLP1R_use2$BMI[which(GLP1R_use2$BMI<bmi1|GLP1R_use2$BMI>bmi2)]<-NA

#T2D
table(GLP1R_use2$T2D,exclude=NULL)

#mood swing
GLP1R_use2$mood_swing[GLP1R_use2$X1920.0.0==0]<-0
GLP1R_use2$mood_swing[GLP1R_use2$X1920.0.0==1]<-1
table(GLP1R_use2$mood_swing,exclude=NULL)

GLP1R_use2$mood_swing2[GLP1R_use2$X1920.0.0==0 & GLP1R_use2$X1940.0.0==0]<-0
GLP1R_use2$mood_swing2[GLP1R_use2$X1920.0.0==1]<-1
table(GLP1R_use2$mood_swing2,exclude=NULL)

#irritable mood
GLP1R_use2$irritable[GLP1R_use2$X1940.0.0==0]<-0
GLP1R_use2$irritable[GLP1R_use2$X1940.0.0==1]<-1
table(GLP1R_use2$irritable,exclude=NULL)

GLP1R_use2$irritable2[GLP1R_use2$X1920.0.0==0 & GLP1R_use2$X1940.0.0==0]<-0
GLP1R_use2$irritable2[GLP1R_use2$X1940.0.0==1]<-1
table(GLP1R_use2$irritable2,exclude=NULL)
################################################################################
#PRS
GLP1R_use2$PRS_trans<-GLP1R_use2$rs10922098*0.101+GLP1R_use2$rs9501393*0.073
GLP1R_use2$PRS_cis_trans<-GLP1R_use2$rs10305432*0.0254+GLP1R_use2$rs10922098*0.101
GLP1R_use2$PRS_3SNPs<-GLP1R_use2$rs10922098*0.101+GLP1R_use2$rs9501393*0.073+GLP1R_use2$rs10305432*0.0254
hist(GLP1R_use2$PRS_trans)
hist(GLP1R_use2$PRS_cis_trans)
hist(GLP1R_use2$PRS_3SNPs)

#SNP-phenotype association
outcomelist<-c("self_harm","suicide","MHS","hospital","GP","help","mood_swing","irritable","mood_swing2","irritable2","T2D")
snplist<-information$rsid
snplist2<-c("PRS_trans","PRS_cis_trans","PRS_3SNPs")
snplist3<-c(snplist,snplist2)

glp1r_binary<-function(pheno,geno){
  eaf <- round(mean(GLP1R_use2[,geno])/2,3)
  #not adjust for age and sex
  firststage<-glm(GLP1R_use2[,pheno] ~., data=GLP1R_use2[,c(geno,"array","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")],binomial(link = "logit"))
  beta1<-summary(firststage)$coefficient[2,"Estimate"]
  se1<-summary(firststage)$coefficient[2,"Std. Error"]
  p1<-summary(firststage)$coefficient[2,"Pr(>|z|)"]
  lower1<-beta1-1.96*se1
  upper1<-beta1+1.96*se1
  OR1<-round(exp(beta1),3)
  ORlower1<-round(exp(lower1),3)
  ORupper1<-round(exp(upper1),3)
  #further adjust for age and sex
  secondstage<-glm(GLP1R_use2[,pheno] ~., data=GLP1R_use2[,c(geno,"array","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")],binomial(link = "logit"))
  beta2<-summary(secondstage)$coefficient[2,"Estimate"]
  se2<-summary(secondstage)$coefficient[2,"Std. Error"]
  p2<-summary(secondstage)$coefficient[2,"Pr(>|z|)"]
  lower2<-beta2-1.96*se2
  upper2<-beta2+1.96*se2
  OR2<-round(exp(beta2),3)
  ORlower2<-round(exp(lower2),3)
  ORupper2<-round(exp(upper2),3)
  write.table(cbind(pheno,geno,eaf,beta1,se1,p1,lower1,upper1,OR1,ORlower1,ORupper1,beta2,se2,p2,lower2,upper2,OR2,ORlower2,ORupper2),file="/path_to_save_results/GLP1R/results.csv", append=TRUE, quote=FALSE, sep=',',row.names=FALSE, col.names=FALSE)
}

for (i in 1:length(snplist3)){
  for (j in 1:length(outcomelist)){
    glp1r_binary(outcomelist[j],snplist3[i])
  }
}

#continuous BMI
glp1r_c<-function(geno){
  eaf <- mean(GLP1R_use2[,geno])/2
  #not adjust for age and sex
  firststage<-lm(GLP1R_use2[,"BMI"] ~., data=GLP1R_use2[,c(geno,"array","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")])
  beta1<-summary(firststage)$coefficient[2,"Estimate"]
  se1<-summary(firststage)$coefficient[2,"Std. Error"]
  p1<-summary(firststage)$coefficient[2,"Pr(>|t|)"]
  lower1<-round(beta1-1.96*se1,3)
  upper1<-round(beta1+1.96*se1,3)
  #further adjust for age and sex
  secondstage<-glm(GLP1R_use2[,"BMI"] ~., data=GLP1R_use2[,c(geno,"array","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")])
  beta2<-summary(secondstage)$coefficient[2,"Estimate"]
  se2<-summary(secondstage)$coefficient[2,"Std. Error"]
  p2<-summary(secondstage)$coefficient[2,"Pr(>|t|)"]
  lower2<-round(beta2-1.96*se2,3)
  upper2<-round(beta2+1.96*se2,3)
  write.table(cbind("BMI",geno,eaf,beta1,se1,p1,lower1,upper1,beta2,se2,p2,lower2,upper2),file="/path_to_save_results/GLP1R/results_BMI.csv", append=TRUE, quote=FALSE, sep=',',row.names=FALSE, col.names=FALSE)
}

for (i in 1:length(snplist3)){glp1r_c(snplist3[i])}

################################################################################
#sex-specific
GLP1R_use2_female<-GLP1R_use2[which(GLP1R_use2$sex==0),]
GLP1R_use2_male<-GLP1R_use2[which(GLP1R_use2$sex==1),]

#SNP-phenotype association
outcomelist<-c("self_harm","suicide","MHS","hospital","GP","help","mood_swing","irritable","mood_swing2","irritable2","T2D")
snplist<-information$rsid
snplist2<-c("PRS_trans","PRS_cis_trans","PRS_3SNPs")
snplist3<-c(snplist,snplist2)
################################################################################
#in female
glp1r_binary_f<-function(pheno,geno){
  eaf <- round(mean(GLP1R_use2_female[,geno])/2,3)
  #not adjust for age and sex
  firststage<-glm(GLP1R_use2_female[,pheno] ~., data=GLP1R_use2_female[,c(geno,"array","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")],binomial(link = "logit"))
  beta1<-summary(firststage)$coefficient[2,"Estimate"]
  se1<-summary(firststage)$coefficient[2,"Std. Error"]
  p1<-summary(firststage)$coefficient[2,"Pr(>|z|)"]
  lower1<-beta1-1.96*se1
  upper1<-beta1+1.96*se1
  OR1<-round(exp(beta1),3)
  ORlower1<-round(exp(lower1),3)
  ORupper1<-round(exp(upper1),3)
  #further adjust for age and sex
  secondstage<-glm(GLP1R_use2_female[,pheno] ~., data=GLP1R_use2_female[,c(geno,"array","age","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")],binomial(link = "logit"))
  beta2<-summary(secondstage)$coefficient[2,"Estimate"]
  se2<-summary(secondstage)$coefficient[2,"Std. Error"]
  p2<-summary(secondstage)$coefficient[2,"Pr(>|z|)"]
  lower2<-beta2-1.96*se2
  upper2<-beta2+1.96*se2
  OR2<-round(exp(beta2),3)
  ORlower2<-round(exp(lower2),3)
  ORupper2<-round(exp(upper2),3)
  write.table(cbind(pheno,geno,eaf,beta1,se1,p1,lower1,upper1,OR1,ORlower1,ORupper1,beta2,se2,p2,lower2,upper2,OR2,ORlower2,ORupper2),file="/path_to_save_results/GLP1R/results_female.csv", append=TRUE, quote=FALSE, sep=',',row.names=FALSE, col.names=FALSE)
}

for (i in 1:length(snplist3)){
  for (j in 1:length(outcomelist)){
    glp1r_binary_f(outcomelist[j],snplist3[i])
  }
}

#in male
glp1r_binary_m<-function(pheno,geno){
  eaf <- round(mean(GLP1R_use2_male[,geno])/2,3)
  #not adjust for age and sex
  firststage<-glm(GLP1R_use2_male[,pheno] ~., data=GLP1R_use2_male[,c(geno,"array","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")],binomial(link = "logit"))
  beta1<-summary(firststage)$coefficient[2,"Estimate"]
  se1<-summary(firststage)$coefficient[2,"Std. Error"]
  p1<-summary(firststage)$coefficient[2,"Pr(>|z|)"]
  lower1<-beta1-1.96*se1
  upper1<-beta1+1.96*se1
  OR1<-round(exp(beta1),3)
  ORlower1<-round(exp(lower1),3)
  ORupper1<-round(exp(upper1),3)
  #further adjust for age and sex
  secondstage<-glm(GLP1R_use2_male[,pheno] ~., data=GLP1R_use2_male[,c(geno,"array","age","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")],binomial(link = "logit"))
  beta2<-summary(secondstage)$coefficient[2,"Estimate"]
  se2<-summary(secondstage)$coefficient[2,"Std. Error"]
  p2<-summary(secondstage)$coefficient[2,"Pr(>|z|)"]
  lower2<-beta2-1.96*se2
  upper2<-beta2+1.96*se2
  OR2<-round(exp(beta2),3)
  ORlower2<-round(exp(lower2),3)
  ORupper2<-round(exp(upper2),3)
  write.table(cbind(pheno,geno,eaf,beta1,se1,p1,lower1,upper1,OR1,ORlower1,ORupper1,beta2,se2,p2,lower2,upper2,OR2,ORlower2,ORupper2),file="/path_to_save_results/GLP1R/results_male.csv", append=TRUE, quote=FALSE, sep=',',row.names=FALSE, col.names=FALSE)
}

for (i in 1:length(snplist3)){
  for (j in 1:length(outcomelist)){
    glp1r_binary_m(outcomelist[j],snplist3[i])
  }
}
################################################################################
#continuous BMI
#in female
glp1r_c_f<-function(geno){
  eaf <- mean(GLP1R_use2_female[,geno])/2
  #not adjust for age and sex
  firststage<-lm(GLP1R_use2_female[,"BMI"] ~., data=GLP1R_use2_female[,c(geno,"array","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")])
  beta1<-summary(firststage)$coefficient[2,"Estimate"]
  se1<-summary(firststage)$coefficient[2,"Std. Error"]
  p1<-summary(firststage)$coefficient[2,"Pr(>|t|)"]
  lower1<-round(beta1-1.96*se1,3)
  upper1<-round(beta1+1.96*se1,3)
  #further adjust for age and sex
  secondstage<-glm(GLP1R_use2_female[,"BMI"] ~., data=GLP1R_use2_female[,c(geno,"array","age","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")])
  beta2<-summary(secondstage)$coefficient[2,"Estimate"]
  se2<-summary(secondstage)$coefficient[2,"Std. Error"]
  p2<-summary(secondstage)$coefficient[2,"Pr(>|t|)"]
  lower2<-round(beta2-1.96*se2,3)
  upper2<-round(beta2+1.96*se2,3)
  write.table(cbind("BMI",geno,eaf,beta1,se1,p1,lower1,upper1,beta2,se2,p2,lower2,upper2),file="/path_to_save_results/GLP1R/results_BMI_female.csv", append=TRUE, quote=FALSE, sep=',',row.names=FALSE, col.names=FALSE)
}

for (i in 1:length(snplist3)){glp1r_c_f(snplist3[i])}

#in male
glp1r_c_m<-function(geno){
  eaf <- mean(GLP1R_use2_male[,geno])/2
  #not adjust for age and sex
  firststage<-lm(GLP1R_use2_male[,"BMI"] ~., data=GLP1R_use2_male[,c(geno,"array","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")])
  beta1<-summary(firststage)$coefficient[2,"Estimate"]
  se1<-summary(firststage)$coefficient[2,"Std. Error"]
  p1<-summary(firststage)$coefficient[2,"Pr(>|t|)"]
  lower1<-round(beta1-1.96*se1,3)
  upper1<-round(beta1+1.96*se1,3)
  #further adjust for age and sex
  secondstage<-glm(GLP1R_use2_male[,"BMI"] ~., data=GLP1R_use2_male[,c(geno,"array","age","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")])
  beta2<-summary(secondstage)$coefficient[2,"Estimate"]
  se2<-summary(secondstage)$coefficient[2,"Std. Error"]
  p2<-summary(secondstage)$coefficient[2,"Pr(>|t|)"]
  lower2<-round(beta2-1.96*se2,3)
  upper2<-round(beta2+1.96*se2,3)
  write.table(cbind("BMI",geno,eaf,beta1,se1,p1,lower1,upper1,beta2,se2,p2,lower2,upper2),file="/path_to_save_results/GLP1R/results_BMI_male.csv", append=TRUE, quote=FALSE, sep=',',row.names=FALSE, col.names=FALSE)
}

for (i in 1:length(snplist3)){glp1r_c_m(snplist3[i])}
