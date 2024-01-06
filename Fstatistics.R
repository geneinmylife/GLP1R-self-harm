###Fstatistics###
exposure_dat<-read.csv("/Users/tq20202/Desktop/ins.csv")
exposure_dat <- cbind(exposure_dat,fstatistics=1)
for (s in 1:nrow(exposure_dat)){
  z <- exposure_dat[s,"beta"]/exposure_dat[s,"se"]
  pve <- z^2/(exposure_dat[s,"samplesize"]+z^2)
  exposure_dat[s,"fstatistics"] <- (exposure_dat[s,"samplesize"]-2)*pve/(1-pve)
}
