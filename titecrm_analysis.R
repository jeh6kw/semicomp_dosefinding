load("results_titecrm_constantgen.Rdata")

#####CONSTANT
titeconstlong <- data.frame(matrix(nrow=18000,ncol=12))
names(titeconstlong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")

titeconstlong$dlttarget <- rep(c(rep(30,3000),rep(20,3000)),3)
titeconstlong$targetdose <- rep(c(rep(1,1000),rep(3,1000),rep(5,1000)),6)
titeconstlong$n <- c(rep(100,6000),rep(60,6000),rep(30,6000))


for(i in 1:18){
  cmat <- resultSet[[i]]
  for(j in 1:1000){
    spec <- cmat[[j]]
    titeconstlong$chosendose[(i-1)*1000+j] <- spec$dose0[nrow(spec)]
    titeconstlong$duration[(i-1)*1000+j] <- spec$accrualtime[nrow(spec)]
    titeconstlong$numdose1[(i-1)*1000+j] <- sum(spec$dose0==0.1)
    titeconstlong$numdose2[(i-1)*1000+j] <- sum(spec$dose0==0.2)
    titeconstlong$numdose3[(i-1)*1000+j] <- sum(spec$dose0==0.4)
    titeconstlong$numdose4[(i-1)*1000+j] <- sum(spec$dose0==0.65)
    titeconstlong$numdose5[(i-1)*1000+j] <- sum(spec$dose0==1)
    titeconstlong$numDLT[(i-1)*1000+j] <- spec$total_dlt[nrow(spec)]
    titeconstlong$numDR[(i-1)*1000+j] <- spec$total_dr[nrow(spec)]
  }
}
titeconstlong$chosendoseorder <- as.integer(as.factor(titeconstlong$chosendose))



titeconstshort <- data.frame(matrix(nrow=18,ncol=16))
names(titeconstshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                           "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
titeconstshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),3)
titeconstshort$targetdose <- rep(c(1,3,5),6)
titeconstshort$n <- c(rep(100,6),rep(60,6),rep(30,6))

durationdiff <- vector()
for(k in 1:18){
  durationdiff[1]<- titeconstlong$duration[((k-1)*1000+1)]
  durationdiff[2:1000] <- diff(titeconstlong$duration[((k-1)*1000+1):(k*1000)],lag=1)
  titeconstshort$avgduration[k] <- mean(durationdiff)
  
  titeconstshort$dose1_selec[k]<-sum(titeconstlong$chosendoseorder[((k-1)*1000+1):(k*1000)]==1)/10
  titeconstshort$dose2_selec[k]<-sum(titeconstlong$chosendoseorder[((k-1)*1000+1):(k*1000)]==2)/10
  titeconstshort$dose3_selec[k]<-sum(titeconstlong$chosendoseorder[((k-1)*1000+1):(k*1000)]==3)/10
  titeconstshort$dose4_selec[k]<-sum(titeconstlong$chosendoseorder[((k-1)*1000+1):(k*1000)]==4)/10
  titeconstshort$dose5_selec[k]<-sum(titeconstlong$chosendoseorder[((k-1)*1000+1):(k*1000)]==5)/10
  
  titeconstshort$dose1_treat[k]<- sum(titeconstlong$numdose1[((k-1)*1000+1):(k*1000)])*100/sum(titeconstlong$n[((k-1)*1000+1):(k*1000)])
  titeconstshort$dose2_treat[k]<- sum(titeconstlong$numdose2[((k-1)*1000+1):(k*1000)])*100/sum(titeconstlong$n[((k-1)*1000+1):(k*1000)])
  titeconstshort$dose3_treat[k]<- sum(titeconstlong$numdose3[((k-1)*1000+1):(k*1000)])*100/sum(titeconstlong$n[((k-1)*1000+1):(k*1000)])
  titeconstshort$dose4_treat[k]<- sum(titeconstlong$numdose4[((k-1)*1000+1):(k*1000)])*100/sum(titeconstlong$n[((k-1)*1000+1):(k*1000)])
  titeconstshort$dose5_treat[k]<- sum(titeconstlong$numdose5[((k-1)*1000+1):(k*1000)])*100/sum(titeconstlong$n[((k-1)*1000+1):(k*1000)])
  
  titeconstshort$avgDLT[k] <- mean(titeconstlong$numDLT[((k-1)*1000+1):(k*1000)])
  titeconstshort$avgDR[k] <- mean(titeconstlong$numDR[((k-1)*1000+1):(k*1000)])
}

write.csv(titeconstshort,"titeconst_summaryresults.csv")




##### WEIBULL
for(file in 1:10){
  
load(paste("results_titecrm_weibullgen_",file,".Rdata",sep=''))


titeweiblong <- data.frame(matrix(nrow=5400,ncol=12))
names(titeweiblong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")

titeweiblong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
titeweiblong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
titeweiblong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))


for(i in 1:54){
  wmat <- resultSet[[i]]
  for(j in 1:100){
    spec <- wmat[[j]]
    titeweiblong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)-1]
    titeweiblong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
    titeweiblong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
    titeweiblong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
    titeweiblong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
    titeweiblong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
    titeweiblong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
    titeweiblong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
    titeweiblong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
  }
}
titeweiblong$chosendoseorder <- as.integer(as.factor(titeweiblong$chosendose))



titeweibshort <- data.frame(matrix(nrow=54,ncol=16))
names(titeweibshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                           "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
titeweibshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
titeweibshort$targetdose <- rep(c(1,3,5),18)
titeweibshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

durationdiff <- vector()
for(k in 1:54){
  durationdiff[1]<- titeweiblong$duration[((k-1)*100+1)]
  durationdiff[2:100] <- diff(titeweiblong$duration[((k-1)*100+1):(k*100)],lag=1)
  titeweibshort$avgduration[k] <- mean(durationdiff)
  
  titeweibshort$dose1_selec[k]<-sum(titeweiblong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
  titeweibshort$dose2_selec[k]<-sum(titeweiblong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
  titeweibshort$dose3_selec[k]<-sum(titeweiblong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
  titeweibshort$dose4_selec[k]<-sum(titeweiblong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
  titeweibshort$dose5_selec[k]<-sum(titeweiblong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
  
  titeweibshort$dose1_treat[k]<- sum(titeweiblong$numdose1[((k-1)*100+1):(k*100)])*100/sum(titeweiblong$n[((k-1)*100+1):(k*100)])
  titeweibshort$dose2_treat[k]<- sum(titeweiblong$numdose2[((k-1)*100+1):(k*100)])*100/sum(titeweiblong$n[((k-1)*100+1):(k*100)])
  titeweibshort$dose3_treat[k]<- sum(titeweiblong$numdose3[((k-1)*100+1):(k*100)])*100/sum(titeweiblong$n[((k-1)*100+1):(k*100)])
  titeweibshort$dose4_treat[k]<- sum(titeweiblong$numdose4[((k-1)*100+1):(k*100)])*100/sum(titeweiblong$n[((k-1)*100+1):(k*100)])
  titeweibshort$dose5_treat[k]<- sum(titeweiblong$numdose5[((k-1)*100+1):(k*100)])*100/sum(titeweiblong$n[((k-1)*100+1):(k*100)])
  titeweibshort$avgDLT[k] <- mean(titeweiblong$numDLT[((k-1)*100+1):(k*100)])
  titeweibshort$avgDR[k] <- mean(titeweiblong$numDR[((k-1)*100+1):(k*100)])
  
}


write.csv(titeweibshort,paste("titeweib_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull
t1 <- read.csv("titeweib_summaryresults_1.csv")
t2<- read.csv("titeweib_summaryresults_2.csv")
t3<- read.csv("titeweib_summaryresults_3.csv")
t4<- read.csv("titeweib_summaryresults_4.csv")
t5<- read.csv("titeweib_summaryresults_5.csv")
t6 <- read.csv("titeweib_summaryresults_6.csv")
t7<- read.csv("titeweib_summaryresults_7.csv")
t8<- read.csv("titeweib_summaryresults_8.csv")
t9<- read.csv("titeweib_summaryresults_9.csv")
t10<- read.csv("titeweib_summaryresults_10.csv")


titeweibshort <- data.frame(matrix(nrow=54,ncol=16))
names(titeweibshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                          "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
titeweibshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
titeweibshort$targetdose <- rep(c(1,3,5),18)
titeweibshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

titeweibshort$avgduration <- (t1$avgduration+t2$avgduration+t3$avgduration+t4$avgduration+t5$avgduration+
                                t6$avgduration+t7$avgduration+t8$avgduration+t9$avgduration+t10$avgduration)/10

titeweibshort$avgDR <- (t1$avgDR+t2$avgDR+t3$avgDR+t4$avgDR+t5$avgDR+t6$avgDR+t7$avgDR+t8$avgDR+t9$avgDR+t10$avgDR)/10

titeweibshort$avgDLT <- (t1$avgDLT+t2$avgDLT+t3$avgDLT+t4$avgDLT+t5$avgDLT+t6$avgDLT+t7$avgDLT+t8$avgDLT+t9$avgDLT+t10$avgDLT)/10

titeweibshort$dose1_selec<-(t1$dose1_selec+t2$dose1_selec+t3$dose1_selec+t4$dose1_selec+t5$dose1_selec+
                              t6$dose1_selec+t7$dose1_selec+t8$dose1_selec+t9$dose1_selec+t10$dose1_selec)/10
titeweibshort$dose2_selec<-(t1$dose2_selec+t2$dose2_selec+t3$dose2_selec+t4$dose2_selec+t5$dose2_selec+
                              t6$dose2_selec+t7$dose2_selec+t8$dose2_selec+t9$dose2_selec+t10$dose2_selec)/10
titeweibshort$dose3_selec<-(t1$dose3_selec+t2$dose3_selec+t3$dose3_selec+t4$dose3_selec+t5$dose3_selec+
                              t6$dose3_selec+t7$dose3_selec+t8$dose3_selec+t9$dose3_selec+t10$dose3_selec)/10
titeweibshort$dose4_selec<-(t1$dose4_selec+t2$dose4_selec+t3$dose4_selec+t4$dose4_selec+t5$dose4_selec+
                              t6$dose4_selec+t7$dose4_selec+t8$dose4_selec+t9$dose4_selec+t10$dose4_selec)/10
titeweibshort$dose5_selec<-(t1$dose5_selec+t2$dose5_selec+t3$dose5_selec+t4$dose5_selec+t5$dose5_selec+
                              t6$dose5_selec+t7$dose5_selec+t8$dose5_selec+t9$dose5_selec+t10$dose5_selec)/10

titeweibshort$dose1_treat<- (t1$dose1_treat+t2$dose1_treat+t3$dose1_treat+t4$dose1_treat+t5$dose1_treat+
                               t6$dose1_treat+t7$dose1_treat+t8$dose1_treat+t9$dose1_treat+t10$dose1_treat)/10
titeweibshort$dose2_treat<- (t1$dose2_treat+t2$dose2_treat+t3$dose2_treat+t4$dose2_treat+t5$dose2_treat+
                               t6$dose2_treat+t7$dose2_treat+t8$dose2_treat+t9$dose2_treat+t10$dose2_treat)/10
titeweibshort$dose3_treat<- (t1$dose3_treat+t2$dose3_treat+t3$dose3_treat+t4$dose3_treat+t5$dose3_treat+
                               t6$dose3_treat+t7$dose3_treat+t8$dose3_treat+t9$dose3_treat+t10$dose3_treat)/10
titeweibshort$dose4_treat<- (t1$dose4_treat+t2$dose4_treat+t3$dose4_treat+t4$dose4_treat+t5$dose4_treat+
                               t6$dose4_treat+t7$dose4_treat+t8$dose4_treat+t9$dose4_treat+t10$dose4_treat)/10
titeweibshort$dose5_treat<- (t1$dose5_treat+t2$dose5_treat+t3$dose5_treat+t4$dose5_treat+t5$dose5_treat+
                               t6$dose5_treat+t7$dose5_treat+t8$dose5_treat+t9$dose5_treat+t10$dose5_treat)/10

write.csv(titeweibshort,file="titeweib_summaryresults.csv")
