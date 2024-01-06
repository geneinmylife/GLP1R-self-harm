#####paiewise zscore test#######
setwd("/Users/tq20202/Desktop")
m1<-read.csv("m1.csv")
m1$diff=m1$beta_F - m1$beta_M
m1$se =sqrt(m1$se_F^2+m1$se_M^2)
m1$Z=m1$diff/m1$se

pvalue.extreme <- function(z) {
  log.pvalue <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
  log10.pvalue <- log.pvalue/log(10) ## from natural log to log10
  mantissa <- 10^(log10.pvalue %% 1)
  exponent <- log10.pvalue %/% 1
  p<-mantissa*10^(exponent)
  return(p)
  ## or return(c(mantissa,exponent))
  #return(sprintf("p value is %1.2f times 10^(%d)",mantissa,exponent))
}

m1$p2<- pvalue.extreme(m1$Z)
write.table(m1, file="ztest_m1.txt", row.names=F, col.names=T, sep="\t", quote=F)

