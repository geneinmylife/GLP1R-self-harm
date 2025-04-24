#package
require(data.table)

#input data
withdraw_list<-read.csv("/path_to_dataset/")
phenotype<-fread("/path_to_dataset/")
glp1r<-fread("/path_to_dataset/")

glp1r_phenotype<-merge(glp1r,phenotype,by="eid")
glp1r_phenotype2<-glp1r_phenotype[-which(glp1r_phenotype$eid%in%withdraw_list$V1),]

#remove non-EUR
table(glp1r_phenotype2$p21000_i0,exclude=NULL)
glp1r_phenotype3<-glp1r_phenotype2[which(glp1r_phenotype2$p21000_i0%in%c("British","Any other white background","Irish","White")),]

#clean outcomes
##ever attempted suicide
glp1r_phenotype3$suicide[glp1r_phenotype3$p20480=="No"]<-0
glp1r_phenotype3$suicide[glp1r_phenotype3$p20483=="Yes"]<-1
table(glp1r_phenotype3$suicide,exclude=NULL)

##ever self-harmed 
glp1r_phenotype3$self_harm[glp1r_phenotype3$p20480=="Yes"]<-1
glp1r_phenotype3$self_harm[glp1r_phenotype3$p20480=="No"]<-0
table(glp1r_phenotype3$self_harm,exclude=NULL)

##mental health services
glp1r_phenotype3$MHS[glp1r_phenotype3$p20554%in%unique(glp1r_phenotype3$p20554)]<-1
glp1r_phenotype3$MHS[glp1r_phenotype3$p20554==""|glp1r_phenotype3$p20554=="Prefer not to answer"]<-NA
glp1r_phenotype3$MHS[is.na(glp1r_phenotype3$MHS) & glp1r_phenotype3$p20480=="No"]<-0
table(glp1r_phenotype3$MHS,exclude=NULL)

##hospital treatment
hospital_treatment<-c("See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|Receive help from friends / family / neighbours",
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|See own GP|Receive help from friends / family / neighbours", 
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|Use a helpline / voluntary organization",                                                           
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)",
"Need hospital treatment (eg A&E)",
"Need hospital treatment (eg A&E)|Receive help from friends / family / neighbours",
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|See own GP",
"Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|See own GP|Receive help from friends / family / neighbours",
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|See own GP|Receive help from friends / family / neighbours",
"Need hospital treatment (eg A&E)|See own GP|Receive help from friends / family / neighbours",                                                                                                       
"Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|Receive help from friends / family / neighbours",                                                                    
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|See own GP",                                                                                  
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|Receive help from friends / family / neighbours",           
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|See own GP",                                             
"Need hospital treatment (eg A&E)|Use a helpline / voluntary organization")

glp1r_phenotype3$hospital[glp1r_phenotype3$p20554%in%hospital_treatment]<-1
glp1r_phenotype3$hospital[glp1r_phenotype3$p20554==""|glp1r_phenotype3$p20554=="Prefer not to answer"]<-NA
glp1r_phenotype3$hospital[glp1r_phenotype3$MHS==0]<-0
table(glp1r_phenotype3$hospital,exclude=NULL)

##see own GP
gp<-c("See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|See own GP|Receive help from friends / family / neighbours",
"See anyone from psychiatric or mental health services, including liaison services|Use a helpline / voluntary organization|See own GP",
"See own GP",                                                                                                                                                                                                       
"See own GP|Receive help from friends / family / neighbours",                                                                                                                                                           
"See anyone from psychiatric or mental health services, including liaison services|See own GP",                                                                                                                       
"Need hospital treatment (eg A&E)|See own GP",                                                                                                                              
"Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|See own GP|Receive help from friends / family / neighbours",                                                                                 
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|See own GP|Receive help from friends / family / neighbours",
"Need hospital treatment (eg A&E)|See own GP|Receive help from friends / family / neighbours",                                                                                                                          
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|See own GP",
"Use a helpline / voluntary organization|See own GP|Receive help from friends / family / neighbours",                                                                                                                 
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|See own GP",
"See anyone from psychiatric or mental health services, including liaison services|Use a helpline / voluntary organization|See own GP|Receive help from friends / family / neighbours",                                
"See anyone from psychiatric or mental health services, including liaison services|See own GP|Receive help from friends / family / neighbours",
"Use a helpline / voluntary organization|See own GP")

glp1r_phenotype3$GP[glp1r_phenotype3$p20554%in%gp]<-1
glp1r_phenotype3$GP[glp1r_phenotype3$p20554==""|glp1r_phenotype3$p20554=="Prefer not to answer"]<-NA
glp1r_phenotype3$GP[glp1r_phenotype3$MHS==0]<-0
table(glp1r_phenotype3$GP,exclude=NULL)

##receive help
receive_help<-c("See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|Receive help from friends / family / neighbours",                                                  
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|See own GP|Receive help from friends / family / neighbours",
"Receive help from friends / family / neighbours",                                                                                                                                                                    
"Need hospital treatment (eg A&E)|Receive help from friends / family / neighbours",                                                                                                                               
"See own GP|Receive help from friends / family / neighbours",                                                                                             
"Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|See own GP|Receive help from friends / family / neighbours",                                                             
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|See own GP|Receive help from friends / family / neighbours",
"Need hospital treatment (eg A&E)|See own GP|Receive help from friends / family / neighbours",                                       
"Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|Receive help from friends / family / neighbours",                                                                                             
"Use a helpline / voluntary organization|Receive help from friends / family / neighbours",                                                                                                                              
"See anyone from psychiatric or mental health services, including liaison services|Use a helpline / voluntary organization|Receive help from friends / family / neighbours",                                         
"See anyone from psychiatric or mental health services, including liaison services|Need hospital treatment (eg A&E)|Use a helpline / voluntary organization|Receive help from friends / family / neighbours",      
"Use a helpline / voluntary organization|See own GP|Receive help from friends / family / neighbours",         
"See anyone from psychiatric or mental health services, including liaison services|Receive help from friends / family / neighbours",                                                                                    
"See anyone from psychiatric or mental health services, including liaison services|Use a helpline / voluntary organization|See own GP|Receive help from friends / family / neighbours",                             
"See anyone from psychiatric or mental health services, including liaison services|See own GP|Receive help from friends / family / neighbours") 
glp1r_phenotype3$help[glp1r_phenotype3$p20554%in%receive_help]<-1
glp1r_phenotype3$help[glp1r_phenotype3$p20554==""|glp1r_phenotype3$p20554=="Prefer not to answer"]<-NA
glp1r_phenotype3$help[glp1r_phenotype3$MHS==0]<-0
table(glp1r_phenotype3$help,exclude=NULL)

##mood swings
glp1r_phenotype3$mood_swing[glp1r_phenotype3$p1920_i0=="No"]<-0
glp1r_phenotype3$mood_swing[glp1r_phenotype3$p1920_i0=="Yes"]<-1
table(glp1r_phenotype3$mood_swing,exclude=NULL)

glp1r_phenotype3$mood_swing2[glp1r_phenotype3$p1920_i0=="No" & glp1r_phenotype3$p1940_i0=="No"]<-0
glp1r_phenotype3$mood_swing2[glp1r_phenotype3$p1920_i0=="Yes"]<-1
table(glp1r_phenotype3$mood_swing2,exclude=NULL)

##irritable mood
glp1r_phenotype3$irritable[glp1r_phenotype3$p1940_i0=="No"]<-0
glp1r_phenotype3$irritable[glp1r_phenotype3$p1940_i0=="Yes"]<-1
table(glp1r_phenotype3$irritable,exclude=NULL)

glp1r_phenotype3$irritable2[glp1r_phenotype3$p1920_i0=="No" & glp1r_phenotype3$p1940_i0=="No"]<-0
glp1r_phenotype3$irritable2[glp1r_phenotype3$p1940_i0=="Yes"]<-1
table(glp1r_phenotype3$irritable2,exclude=NULL)

#clean covariates
##model 1
###age
summary(glp1r_phenotype3$p21022)

###sex(0=female, 1=male)
table(glp1r_phenotype3$p31,exclude=NULL)
glp1r_phenotype3$sex[glp1r_phenotype3$p31=="Female"]<-0
glp1r_phenotype3$sex[glp1r_phenotype3$p31=="Male"]<-1
table(glp1r_phenotype3$sex,exclude=NULL)

###UKB assessment center
table(glp1r_phenotype3$p54_i0,exclude=NULL)

###Townsend deprivation index
summary(glp1r_phenotype3$p22189)

###education
unique(glp1r_phenotype3$p6138_i0)
university<-c("College or University degree",
"College or University degree|A levels/AS levels or equivalent|O levels/GCSEs or equivalent|Other professional qualifications eg: nursing, teaching",
"College or University degree|Other professional qualifications eg: nursing, teaching",
"College or University degree|A levels/AS levels or equivalent|O levels/GCSEs or equivalent|CSEs or equivalent|Other professional qualifications eg: nursing, teaching",
"College or University degree|A levels/AS levels or equivalent|O levels/GCSEs or equivalent",                                                                                                
"College or University degree|NVQ or HND or HNC or equivalent|Other professional qualifications eg: nursing, teaching",                                                                                 
"College or University degree|O levels/GCSEs or equivalent|NVQ or HND or HNC or equivalent",                                                                                                            
"College or University degree|O levels/GCSEs or equivalent",                                                                                                                                            
"College or University degree|A levels/AS levels or equivalent|O levels/GCSEs or equivalent|NVQ or HND or HNC or equivalent|Other professional qualifications eg: nursing, teaching",                   
"College or University degree|A levels/AS levels or equivalent|O levels/GCSEs or equivalent|CSEs or equivalent|NVQ or HND or HNC or equivalent",                                                  
"College or University degree|A levels/AS levels or equivalent|O levels/GCSEs or equivalent|CSEs or equivalent",                                                                                        
"College or University degree|O levels/GCSEs or equivalent|CSEs or equivalent|Other professional qualifications eg: nursing, teaching",                                                                 
"College or University degree|A levels/AS levels or equivalent|Other professional qualifications eg: nursing, teaching",                                                                                
"College or University degree|A levels/AS levels or equivalent",                                                                                                                                        
"College or University degree|O levels/GCSEs or equivalent|Other professional qualifications eg: nursing, teaching",                                                                                    
"College or University degree|A levels/AS levels or equivalent|O levels/GCSEs or equivalent|NVQ or HND or HNC or equivalent",                                                                           
"College or University degree|A levels/AS levels or equivalent|O levels/GCSEs or equivalent|CSEs or equivalent|NVQ or HND or HNC or equivalent|Other professional qualifications eg: nursing, teaching",
"College or University degree|O levels/GCSEs or equivalent|CSEs or equivalent|NVQ or HND or HNC or equivalent|Other professional qualifications eg: nursing, teaching",
"College or University degree|O levels/GCSEs or equivalent|CSEs or equivalent|NVQ or HND or HNC or equivalent",
"College or University degree|NVQ or HND or HNC or equivalent",
"College or University degree|CSEs or equivalent|Other professional qualifications eg: nursing, teaching",                                                                                              
"College or University degree|CSEs or equivalent",
"College or University degree|A levels/AS levels or equivalent|NVQ or HND or HNC or equivalent|Other professional qualifications eg: nursing, teaching",                                                
"College or University degree|O levels/GCSEs or equivalent|NVQ or HND or HNC or equivalent|Other professional qualifications eg: nursing, teaching",
"College or University degree|A levels/AS levels or equivalent|NVQ or HND or HNC or equivalent",                                                                                                        
"College or University degree|CSEs or equivalent|NVQ or HND or HNC or equivalent",
"College or University degree|O levels/GCSEs or equivalent|CSEs or equivalent",                                                                                                                        
"College or University degree|A levels/AS levels or equivalent|CSEs or equivalent",                                                                                                                     
"College or University degree|A levels/AS levels or equivalent|CSEs or equivalent|Other professional qualifications eg: nursing, teaching",                                                             
"College or University degree|A levels/AS levels or equivalent|CSEs or equivalent|NVQ or HND or HNC or equivalent|Other professional qualifications eg: nursing, teaching",                             
"College or University degree|A levels/AS levels or equivalent|CSEs or equivalent|NVQ or HND or HNC or equivalent",                                                                                     
"College or University degree|CSEs or equivalent|NVQ or HND or HNC or equivalent|Other professional qualifications eg: nursing, teaching")

glp1r_phenotype3$edu<-0
glp1r_phenotype3$edu[glp1r_phenotype3$p6138_i0%in%university]<-1
glp1r_phenotype3$edu[glp1r_phenotype3$p6138_i0%in%c("","None of the above","Prefer not to answer")]<-NA
table(glp1r_phenotype3$edu,exclude=NULL)

###annual household income
table(glp1r_phenotype3$p738_i0,exclude=NULL)
glp1r_phenotype3$income[glp1r_phenotype3$p738_i0=="Less than 18,000"]<-1
glp1r_phenotype3$income[glp1r_phenotype3$p738_i0=="18,000 to 30,999"]<-2
glp1r_phenotype3$income[glp1r_phenotype3$p738_i0=="31,000 to 51,999"]<-3
glp1r_phenotype3$income[glp1r_phenotype3$p738_i0=="52,000 to 100,000"]<-4
glp1r_phenotype3$income[glp1r_phenotype3$p738_i0=="Greater than 100,000"]<-5
table(glp1r_phenotype3$income,exclude=NULL)

###BMI 
summary(glp1r_phenotype3$p21001_i0)

###smoking status
table(glp1r_phenotype3$p20116_i0,exclude=NULL)
glp1r_phenotype3$smoke[glp1r_phenotype3$p20116_i0=="Never"]<-0
glp1r_phenotype3$smoke[glp1r_phenotype3$p20116_i0=="Previous"]<-1
glp1r_phenotype3$smoke[glp1r_phenotype3$p20116_i0=="Current"]<-2
table(glp1r_phenotype3$smoke,exclude=NULL)

###alcohol consumption
table(glp1r_phenotype3$p1558_i0,exclude=NULL)
glp1r_phenotype3$alcohol[glp1r_phenotype3$p1558_i0=="Never"]<-1
glp1r_phenotype3$alcohol[glp1r_phenotype3$p1558_i0=="Special occasions only"]<-2
glp1r_phenotype3$alcohol[glp1r_phenotype3$p1558_i0=="One to three times a month"]<-3
glp1r_phenotype3$alcohol[glp1r_phenotype3$p1558_i0=="Once or twice a week"]<-4
glp1r_phenotype3$alcohol[glp1r_phenotype3$p1558_i0=="Three or four times a week"]<-5
glp1r_phenotype3$alcohol[glp1r_phenotype3$p1558_i0=="Daily or almost daily"]<-6
table(glp1r_phenotype3$alcohol,exclude=NULL)

################################################################################
#logistic regression
##model 1
fit1<-glm(glp1r_phenotype3$suicide ~ glp1r+p21022+sex+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit2<-glm(glp1r_phenotype3$self_harm ~ glp1r+p21022+sex+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit3<-glm(glp1r_phenotype3$MHS ~ glp1r+p21022+sex+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit4<-glm(glp1r_phenotype3$hospital ~ glp1r+p21022+sex+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit5<-glm(glp1r_phenotype3$GP ~ glp1r+p21022+sex+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit6<-glm(glp1r_phenotype3$help ~ glp1r+p21022+sex+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit7<-glm(glp1r_phenotype3$mood_swing ~ glp1r+p21022+sex+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit8<-glm(glp1r_phenotype3$mood_swing2 ~ glp1r+p21022+sex+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit9<-glm(glp1r_phenotype3$irritable ~ glp1r+p21022+sex+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit10<-glm(glp1r_phenotype3$irritable2 ~ glp1r+p21022+sex+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())

extract_results1<-function(model1){
  sumx1<-summary(model1)
  logor1=sumx1$coefficient[2,"Estimate"]
  se1=sumx1$coefficient[2,"Std. Error"]
  p1=sumx1$coefficient[2,"Pr(>|z|)"]
  or1<-round(exp(logor1),3)
  lci1<-round(exp(logor1-1.96*se1),3)
  uci1<-round(exp(logor1+1.96*se1),3)
  write.table(cbind("model1",logor1,se1,p1,or1,lci1,uci1),
              file="/path_to_results/",
              append = TRUE,quote=FALSE,sep=",",col.names = FALSE,row.names =TRUE)
}

extract_results1(fit1)
extract_results1(fit2)
extract_results1(fit3)
extract_results1(fit4)
extract_results1(fit5)
extract_results1(fit6)
extract_results1(fit7)
extract_results1(fit8)
extract_results1(fit9)
extract_results1(fit10)

glp1r_phenotype3_model1<-glp1r_phenotype3[which(!is.na(glp1r_phenotype3$p21022) & 
                                                  !is.na(glp1r_phenotype3$sex) & 
                                                  !is.na(glp1r_phenotype3$p21001_i0) &
                                                  !is.na(glp1r_phenotype3$glp1r)),]
table(glp1r_phenotype3_model1$suicide,exclude=NULL)
table(glp1r_phenotype3_model1$self_harm,exclude=NULL)
table(glp1r_phenotype3_model1$MHS,exclude=NULL)
table(glp1r_phenotype3_model1$hospital,exclude=NULL)
table(glp1r_phenotype3_model1$GP,exclude=NULL)
table(glp1r_phenotype3_model1$help,exclude=NULL)
table(glp1r_phenotype3_model1$mood_swing,exclude=NULL)
table(glp1r_phenotype3_model1$mood_swing2,exclude=NULL)
table(glp1r_phenotype3_model1$irritable,exclude=NULL)
table(glp1r_phenotype3_model1$irritable2,exclude=NULL)

glp1r_phenotype3_model2<-glp1r_phenotype3[which(!is.na(glp1r_phenotype3$p21022) & 
                                                  !is.na(glp1r_phenotype3$sex) & 
                                                  !is.na(glp1r_phenotype3$p21001_i0) &
                                                  !is.na(glp1r_phenotype3$p54_i0) &
                                                  !is.na(glp1r_phenotype3$p22189) & 
                                                  !is.na(glp1r_phenotype3$smoke) &
                                                  !is.na(glp1r_phenotype3$alcohol) &
                                                  !is.na(glp1r_phenotype3$glp1r)),]
table(glp1r_phenotype3_model2$suicide,exclude=NULL)
table(glp1r_phenotype3_model2$self_harm,exclude=NULL)
table(glp1r_phenotype3_model2$MHS,exclude=NULL)
table(glp1r_phenotype3_model2$hospital,exclude=NULL)
table(glp1r_phenotype3_model2$GP,exclude=NULL)
table(glp1r_phenotype3_model2$help,exclude=NULL)
table(glp1r_phenotype3_model2$mood_swing,exclude=NULL)
table(glp1r_phenotype3_model2$mood_swing2,exclude=NULL)
table(glp1r_phenotype3_model2$irritable,exclude=NULL)
table(glp1r_phenotype3_model2$irritable2,exclude=NULL)

##model 2
fit1<-glm(glp1r_phenotype3$suicide ~ glp1r+p21022+sex+p54_i0+p22189+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit2<-glm(glp1r_phenotype3$self_harm ~ glp1r+p21022+sex+p54_i0+p22189+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit3<-glm(glp1r_phenotype3$MHS ~ glp1r+p21022+sex+p54_i0+p22189+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit4<-glm(glp1r_phenotype3$hospital ~ glp1r+p21022+sex+p54_i0+p22189+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit5<-glm(glp1r_phenotype3$GP ~ glp1r+p21022+sex+p54_i0+p22189+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit6<-glm(glp1r_phenotype3$help ~ glp1r+p21022+sex+p54_i0+p22189+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit7<-glm(glp1r_phenotype3$mood_swing ~ glp1r+p21022+sex+p54_i0+p22189+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit8<-glm(glp1r_phenotype3$mood_swing2 ~ glp1r+p21022+sex+p54_i0+p22189+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit9<-glm(glp1r_phenotype3$irritable ~ glp1r+p21022+sex+p54_i0+p22189+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())
fit10<-glm(glp1r_phenotype3$irritable2 ~ glp1r+p21022+sex+p54_i0+p22189+p21001_i0+smoke+alcohol,data=glp1r_phenotype3,family=binomial())

extract_results2<-function(model2){
  sumx1<-summary(model2)
  logor1=sumx1$coefficient[2,"Estimate"]
  se1=sumx1$coefficient[2,"Std. Error"]
  p1=sumx1$coefficient[2,"Pr(>|z|)"]
  or1<-round(exp(logor1),3)
  lci1<-round(exp(logor1-1.96*se1),3)
  uci1<-round(exp(logor1+1.96*se1),3) 
  write.table(cbind("model2",logor1,se1,p1,or1,lci1,uci1),
              file="/path_to_results/",
              append = TRUE,quote=FALSE,sep=",",col.names = FALSE,row.names =TRUE)
}

extract_results2(fit1)
extract_results2(fit2)
extract_results2(fit3)
extract_results2(fit4)
extract_results2(fit5)
extract_results2(fit6)
extract_results2(fit7)
extract_results2(fit8)
extract_results2(fit9)
extract_results2(fit10)

