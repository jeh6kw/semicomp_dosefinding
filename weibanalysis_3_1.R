####################All Weibull 3:1 rule analysis. C, BD, SD, INC gen. All uninformative priors


###Const gen

for(file in c(3,5,6:11)){
  load(paste("results_cgen_weibsis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  cgen_weib31long <- data.frame(matrix(nrow=5400,ncol=12))
  names(cgen_weib31long) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  cgen_weib31long$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  cgen_weib31long$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  cgen_weib31long$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      cgen_weib31long$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      cgen_weib31long$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      cgen_weib31long$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      cgen_weib31long$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      cgen_weib31long$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      cgen_weib31long$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      cgen_weib31long$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      cgen_weib31long$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      cgen_weib31long$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  cgen_weib31long$chosendoseorder <- as.integer(as.factor(cgen_weib31long$chosendose))
  
  
  
  cgen_weib31short <- data.frame(matrix(nrow=54,ncol=16))
  names(cgen_weib31short) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                              "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  cgen_weib31short$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  cgen_weib31short$targetdose <- rep(c(1,3,5),18)
  cgen_weib31short$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- cgen_weib31long$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(cgen_weib31long$duration[((k-1)*100+1):(k*100)],lag=1)
    cgen_weib31short$avgduration[k] <- mean(durationdiff)
    
    cgen_weib31short$dose1_selec[k]<-sum(cgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    cgen_weib31short$dose2_selec[k]<-sum(cgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    cgen_weib31short$dose3_selec[k]<-sum(cgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    cgen_weib31short$dose4_selec[k]<-sum(cgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    cgen_weib31short$dose5_selec[k]<-sum(cgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    cgen_weib31short$dose1_treat[k]<- sum(cgen_weib31long$numdose1[((k-1)*100+1):(k*100)])*100/sum(cgen_weib31long$n[((k-1)*100+1):(k*100)])
    cgen_weib31short$dose2_treat[k]<- sum(cgen_weib31long$numdose2[((k-1)*100+1):(k*100)])*100/sum(cgen_weib31long$n[((k-1)*100+1):(k*100)])
    cgen_weib31short$dose3_treat[k]<- sum(cgen_weib31long$numdose3[((k-1)*100+1):(k*100)])*100/sum(cgen_weib31long$n[((k-1)*100+1):(k*100)])
    cgen_weib31short$dose4_treat[k]<- sum(cgen_weib31long$numdose4[((k-1)*100+1):(k*100)])*100/sum(cgen_weib31long$n[((k-1)*100+1):(k*100)])
    cgen_weib31short$dose5_treat[k]<- sum(cgen_weib31long$numdose5[((k-1)*100+1):(k*100)])*100/sum(cgen_weib31long$n[((k-1)*100+1):(k*100)])
    
    cgen_weib31short$avgDLT[k] <- mean(cgen_weib31long$numDLT[((k-1)*100+1):(k*100)])
    cgen_weib31short$avgDR[k] <- mean(cgen_weib31long$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(cgen_weib31short,paste("cgen_weib31_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("cgen_weib31_summaryresults_11.csv")
c2<- read.csv("cgen_weib31_summaryresults_10.csv")
c3<- read.csv("cgen_weib31_summaryresults_3.csv")
c4<- read.csv("cgen_weib31_summaryresults_9.csv")
c5<- read.csv("cgen_weib31_summaryresults_5.csv")
c6 <- read.csv("cgen_weib31_summaryresults_6.csv")
c7<- read.csv("cgen_weib31_summaryresults_7.csv")
c8<- read.csv("cgen_weib31_summaryresults_8.csv")



cgen_weib31short <- data.frame(matrix(nrow=54,ncol=16))
names(cgen_weib31short) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                            "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
cgen_weib31short$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
cgen_weib31short$targetdose <- rep(c(1,3,5),18)
cgen_weib31short$n <- c(rep(100,18),rep(60,18),rep(30,18))

cgen_weib31short$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                  c6$avgduration+c7$avgduration+c8$avgduration)*100/800

cgen_weib31short$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                             c6$avgDLT+c7$avgDLT+c8$avgDLT)*100/800

cgen_weib31short$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                            c6$avgDR+c7$avgDR+c8$avgDR)*100/800

cgen_weib31short$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                                c6$dose1_selec+c7$dose1_selec+c8$dose1_selec)*100/800
cgen_weib31short$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                                c6$dose2_selec+c7$dose2_selec+c8$dose2_selec)*100/800
cgen_weib31short$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                                c6$dose3_selec+c7$dose3_selec+c8$dose3_selec)*100/800
cgen_weib31short$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                                c6$dose4_selec+c7$dose4_selec+c8$dose4_selec)*100/800
cgen_weib31short$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                                c6$dose5_selec+c7$dose5_selec+c8$dose5_selec)*100/800

cgen_weib31short$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                 c6$dose1_treat+c7$dose1_treat+c8$dose1_treat)*100/800
cgen_weib31short$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                 c6$dose2_treat+c7$dose2_treat+c8$dose2_treat)*100/800
cgen_weib31short$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                 c6$dose3_treat+c7$dose3_treat+c8$dose3_treat)*100/800
cgen_weib31short$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                 c6$dose4_treat+c7$dose4_treat+c8$dose4_treat)*100/800
cgen_weib31short$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                 c6$dose5_treat+c7$dose5_treat+c8$dose5_treat)*100/800



write.csv(cgen_weib31short,file="cgen_weib31_summaryresults.csv")








########################################################################################################
########################################################################################################






### Big Decreasing Gen


for(file in c(1:4,6:11)){
  load(paste("results_wbdgen_weibsis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  bdgen_weib31long <- data.frame(matrix(nrow=5400,ncol=12))
  names(bdgen_weib31long) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  bdgen_weib31long$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  bdgen_weib31long$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  bdgen_weib31long$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      bdgen_weib31long$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      bdgen_weib31long$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      bdgen_weib31long$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      bdgen_weib31long$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      bdgen_weib31long$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      bdgen_weib31long$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      bdgen_weib31long$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      bdgen_weib31long$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      bdgen_weib31long$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  bdgen_weib31long$chosendoseorder <- as.integer(as.factor(bdgen_weib31long$chosendose))
  
  
  
  bdgen_weib31short <- data.frame(matrix(nrow=54,ncol=16))
  names(bdgen_weib31short) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                           "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  bdgen_weib31short$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  bdgen_weib31short$targetdose <- rep(c(1,3,5),18)
  bdgen_weib31short$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- bdgen_weib31long$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(bdgen_weib31long$duration[((k-1)*100+1):(k*100)],lag=1)
    bdgen_weib31short$avgduration[k] <- mean(durationdiff)
    
    bdgen_weib31short$dose1_selec[k]<-sum(bdgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    bdgen_weib31short$dose2_selec[k]<-sum(bdgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    bdgen_weib31short$dose3_selec[k]<-sum(bdgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    bdgen_weib31short$dose4_selec[k]<-sum(bdgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    bdgen_weib31short$dose5_selec[k]<-sum(bdgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    bdgen_weib31short$dose1_treat[k]<- sum(bdgen_weib31long$numdose1[((k-1)*100+1):(k*100)])*100/sum(bdgen_weib31long$n[((k-1)*100+1):(k*100)])
    bdgen_weib31short$dose2_treat[k]<- sum(bdgen_weib31long$numdose2[((k-1)*100+1):(k*100)])*100/sum(bdgen_weib31long$n[((k-1)*100+1):(k*100)])
    bdgen_weib31short$dose3_treat[k]<- sum(bdgen_weib31long$numdose3[((k-1)*100+1):(k*100)])*100/sum(bdgen_weib31long$n[((k-1)*100+1):(k*100)])
    bdgen_weib31short$dose4_treat[k]<- sum(bdgen_weib31long$numdose4[((k-1)*100+1):(k*100)])*100/sum(bdgen_weib31long$n[((k-1)*100+1):(k*100)])
    bdgen_weib31short$dose5_treat[k]<- sum(bdgen_weib31long$numdose5[((k-1)*100+1):(k*100)])*100/sum(bdgen_weib31long$n[((k-1)*100+1):(k*100)])
    
    bdgen_weib31short$avgDLT[k] <- mean(bdgen_weib31long$numDLT[((k-1)*100+1):(k*100)])
    bdgen_weib31short$avgDR[k] <- mean(bdgen_weib31long$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(bdgen_weib31short,paste("bdgen_weib31_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("bdgen_weib31_summaryresults_1.csv")
c2<- read.csv("bdgen_weib31_summaryresults_2.csv")
c3<- read.csv("bdgen_weib31_summaryresults_3.csv")
c4<- read.csv("bdgen_weib31_summaryresults_4.csv")
c5<- read.csv("bdgen_weib31_summaryresults_11.csv")
c6 <- read.csv("bdgen_weib31_summaryresults_6.csv")
c7<- read.csv("bdgen_weib31_summaryresults_7.csv")
c8<- read.csv("bdgen_weib31_summaryresults_8.csv")
c9<- read.csv("bdgen_weib31_summaryresults_9.csv")
c10<- read.csv("bdgen_weib31_summaryresults_10.csv")


bdgen_weib31short <- data.frame(matrix(nrow=54,ncol=16))
names(bdgen_weib31short) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                         "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
bdgen_weib31short$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
bdgen_weib31short$targetdose <- rep(c(1,3,5),18)
bdgen_weib31short$n <- c(rep(100,18),rep(60,18),rep(30,18))

bdgen_weib31short$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                               c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

bdgen_weib31short$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                          c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

bdgen_weib31short$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                         c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

bdgen_weib31short$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                             c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
bdgen_weib31short$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                             c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
bdgen_weib31short$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                             c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
bdgen_weib31short$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                             c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
bdgen_weib31short$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                             c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

bdgen_weib31short$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                              c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
bdgen_weib31short$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                              c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
bdgen_weib31short$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                              c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
bdgen_weib31short$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                              c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
bdgen_weib31short$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                              c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10



write.csv(bdgen_weib31short,file="bdgen_weib31_summaryresults.csv")








########################################################################################################
########################################################################################################




### Small Decreasing Gen

for(file in c(2,3,5,7:9,11,12)){
  load(paste("results_wsdgen_weibsis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  sdgen_weib31long <- data.frame(matrix(nrow=5400,ncol=12))
  names(sdgen_weib31long) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  sdgen_weib31long$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  sdgen_weib31long$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  sdgen_weib31long$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      sdgen_weib31long$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      sdgen_weib31long$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      sdgen_weib31long$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      sdgen_weib31long$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      sdgen_weib31long$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      sdgen_weib31long$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      sdgen_weib31long$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      sdgen_weib31long$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      sdgen_weib31long$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  sdgen_weib31long$chosendoseorder <- as.integer(as.factor(sdgen_weib31long$chosendose))
  
  
  
  sdgen_weib31short <- data.frame(matrix(nrow=54,ncol=16))
  names(sdgen_weib31short) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                           "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  sdgen_weib31short$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  sdgen_weib31short$targetdose <- rep(c(1,3,5),18)
  sdgen_weib31short$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- sdgen_weib31long$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(sdgen_weib31long$duration[((k-1)*100+1):(k*100)],lag=1)
    sdgen_weib31short$avgduration[k] <- mean(durationdiff)
    
    sdgen_weib31short$dose1_selec[k]<-sum(sdgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    sdgen_weib31short$dose2_selec[k]<-sum(sdgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    sdgen_weib31short$dose3_selec[k]<-sum(sdgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    sdgen_weib31short$dose4_selec[k]<-sum(sdgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    sdgen_weib31short$dose5_selec[k]<-sum(sdgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    sdgen_weib31short$dose1_treat[k]<- sum(sdgen_weib31long$numdose1[((k-1)*100+1):(k*100)])*100/sum(sdgen_weib31long$n[((k-1)*100+1):(k*100)])
    sdgen_weib31short$dose2_treat[k]<- sum(sdgen_weib31long$numdose2[((k-1)*100+1):(k*100)])*100/sum(sdgen_weib31long$n[((k-1)*100+1):(k*100)])
    sdgen_weib31short$dose3_treat[k]<- sum(sdgen_weib31long$numdose3[((k-1)*100+1):(k*100)])*100/sum(sdgen_weib31long$n[((k-1)*100+1):(k*100)])
    sdgen_weib31short$dose4_treat[k]<- sum(sdgen_weib31long$numdose4[((k-1)*100+1):(k*100)])*100/sum(sdgen_weib31long$n[((k-1)*100+1):(k*100)])
    sdgen_weib31short$dose5_treat[k]<- sum(sdgen_weib31long$numdose5[((k-1)*100+1):(k*100)])*100/sum(sdgen_weib31long$n[((k-1)*100+1):(k*100)])
    
    sdgen_weib31short$avgDLT[k] <- mean(sdgen_weib31long$numDLT[((k-1)*100+1):(k*100)])
    sdgen_weib31short$avgDR[k] <- mean(sdgen_weib31long$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(sdgen_weib31short,paste("sdgen_weib31_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c11 <- read.csv("sdgen_weib31_summaryresults_11.csv")
c2<- read.csv("sdgen_weib31_summaryresults_2.csv")
c3<- read.csv("sdgen_weib31_summaryresults_3.csv")
#c4<- read.csv("sdgen_weib31_summaryresults_4.csv")
c5<- read.csv("sdgen_weib31_summaryresults_5.csv")
#c6 <- read.csv("sdgen_weib31_summaryresults_6.csv")
c7<- read.csv("sdgen_weib31_summaryresults_7.csv")
c8<- read.csv("sdgen_weib31_summaryresults_8.csv")
c9<- read.csv("sdgen_weib31_summaryresults_9.csv")
c12<- read.csv("sdgen_weib31_summaryresults_12.csv")


sdgen_weib31short <- data.frame(matrix(nrow=54,ncol=16))
names(sdgen_weib31short) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                         "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
sdgen_weib31short$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
sdgen_weib31short$targetdose <- rep(c(1,3,5),18)
sdgen_weib31short$n <- c(rep(100,18),rep(60,18),rep(30,18))

sdgen_weib31short$avgduration <- (c2$avgduration+c3$avgduration+c11$avgduration+c5$avgduration+
                               c12$avgduration+c7$avgduration+c8$avgduration+c9$avgduration)*100/800

sdgen_weib31short$avgDLT <- (c2$avgDLT+c3$avgDLT+c11$avgDLT+c5$avgDLT+
                          c12$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT)*100/800

sdgen_weib31short$avgDR <- (c2$avgDR+c3$avgDR+c11$avgDR+c5$avgDR+
                         c12$avgDR+c7$avgDR+c8$avgDR+c9$avgDR)*100/800

sdgen_weib31short$dose1_selec<-(c2$dose1_selec+c3$dose1_selec+c11$dose1_selec+c5$dose1_selec+
                             c12$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec)*100/800
sdgen_weib31short$dose2_selec<-(c2$dose2_selec+c3$dose2_selec+c11$dose2_selec+c5$dose2_selec+
                             c12$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec)*100/800
sdgen_weib31short$dose3_selec<-(c2$dose3_selec+c3$dose3_selec+c11$dose3_selec+c5$dose3_selec+
                             c12$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec)*100/800
sdgen_weib31short$dose4_selec<-(c2$dose4_selec+c3$dose4_selec+c11$dose4_selec+c5$dose4_selec+
                             c12$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec)*100/800
sdgen_weib31short$dose5_selec<-(c2$dose5_selec+c3$dose5_selec+c11$dose5_selec+c5$dose5_selec+
                             c12$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec)*100/800

sdgen_weib31short$dose1_treat<- (c2$dose1_treat+c3$dose1_treat+c11$dose1_treat+c5$dose1_treat+
                              c12$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat)*100/800
sdgen_weib31short$dose2_treat<- (c2$dose2_treat+c3$dose2_treat+c11$dose2_treat+c5$dose2_treat+
                              c12$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat)*100/800
sdgen_weib31short$dose3_treat<- (c2$dose3_treat+c3$dose3_treat+c11$dose3_treat+c5$dose3_treat+
                              c12$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat)*100/800
sdgen_weib31short$dose4_treat<- (c2$dose4_treat+c3$dose4_treat+c11$dose4_treat+c5$dose4_treat+
                              c12$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat)*100/800
sdgen_weib31short$dose5_treat<- (c2$dose5_treat+c3$dose5_treat+c11$dose5_treat+c5$dose5_treat+
                              c12$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat)*100/800



write.csv(sdgen_weib31short,file="sdgen_weib31_summaryresults.csv")







########################################################################################################
########################################################################################################




#### Increasing gen 



for(file in c(1,2,4:7,10,12)){
  load(paste("results_wincgen_weibsis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  incgen_weib31long <- data.frame(matrix(nrow=5400,ncol=12))
  names(incgen_weib31long) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  incgen_weib31long$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  incgen_weib31long$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  incgen_weib31long$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      incgen_weib31long$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      incgen_weib31long$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      incgen_weib31long$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      incgen_weib31long$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      incgen_weib31long$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      incgen_weib31long$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      incgen_weib31long$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      incgen_weib31long$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      incgen_weib31long$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  incgen_weib31long$chosendoseorder <- as.integer(as.factor(incgen_weib31long$chosendose))
  
  
  
  incgen_weib31short <- data.frame(matrix(nrow=54,ncol=16))
  names(incgen_weib31short) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                            "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  incgen_weib31short$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  incgen_weib31short$targetdose <- rep(c(1,3,5),18)
  incgen_weib31short$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- incgen_weib31long$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(incgen_weib31long$duration[((k-1)*100+1):(k*100)],lag=1)
    incgen_weib31short$avgduration[k] <- mean(durationdiff)
    
    incgen_weib31short$dose1_selec[k]<-sum(incgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    incgen_weib31short$dose2_selec[k]<-sum(incgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    incgen_weib31short$dose3_selec[k]<-sum(incgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    incgen_weib31short$dose4_selec[k]<-sum(incgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    incgen_weib31short$dose5_selec[k]<-sum(incgen_weib31long$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    incgen_weib31short$dose1_treat[k]<- sum(incgen_weib31long$numdose1[((k-1)*100+1):(k*100)])*100/sum(incgen_weib31long$n[((k-1)*100+1):(k*100)])
    incgen_weib31short$dose2_treat[k]<- sum(incgen_weib31long$numdose2[((k-1)*100+1):(k*100)])*100/sum(incgen_weib31long$n[((k-1)*100+1):(k*100)])
    incgen_weib31short$dose3_treat[k]<- sum(incgen_weib31long$numdose3[((k-1)*100+1):(k*100)])*100/sum(incgen_weib31long$n[((k-1)*100+1):(k*100)])
    incgen_weib31short$dose4_treat[k]<- sum(incgen_weib31long$numdose4[((k-1)*100+1):(k*100)])*100/sum(incgen_weib31long$n[((k-1)*100+1):(k*100)])
    incgen_weib31short$dose5_treat[k]<- sum(incgen_weib31long$numdose5[((k-1)*100+1):(k*100)])*100/sum(incgen_weib31long$n[((k-1)*100+1):(k*100)])
    
    incgen_weib31short$avgDLT[k] <- mean(incgen_weib31long$numDLT[((k-1)*100+1):(k*100)])
    incgen_weib31short$avgDR[k] <- mean(incgen_weib31long$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(incgen_weib31short,paste("incgen_weib31_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("incgen_weib31_summaryresults_1.csv")
c2<- read.csv("incgen_weib31_summaryresults_2.csv")
c3<- read.csv("incgen_weib31_summaryresults_10.csv")
c4<- read.csv("incgen_weib31_summaryresults_4.csv")
c5<- read.csv("incgen_weib31_summaryresults_5.csv")
c6 <- read.csv("incgen_weib31_summaryresults_6.csv")
c7<- read.csv("incgen_weib31_summaryresults_7.csv")
c8<- read.csv("incgen_weib31_summaryresults_12.csv")



incgen_weib31short <- data.frame(matrix(nrow=54,ncol=16))
names(incgen_weib31short) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                          "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
incgen_weib31short$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
incgen_weib31short$targetdose <- rep(c(1,3,5),18)
incgen_weib31short$n <- c(rep(100,18),rep(60,18),rep(30,18))

incgen_weib31short$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                c6$avgduration+c7$avgduration+c8$avgduration)*100/800

incgen_weib31short$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                           c6$avgDLT+c7$avgDLT+c8$avgDLT)*100/800

incgen_weib31short$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                          c6$avgDR+c7$avgDR+c8$avgDR)*100/800

incgen_weib31short$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                              c6$dose1_selec+c7$dose1_selec+c8$dose1_selec)*100/800
incgen_weib31short$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                              c6$dose2_selec+c7$dose2_selec+c8$dose2_selec)*100/800
incgen_weib31short$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                              c6$dose3_selec+c7$dose3_selec+c8$dose3_selec)*100/800
incgen_weib31short$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                              c6$dose4_selec+c7$dose4_selec+c8$dose4_selec)*100/800
incgen_weib31short$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                              c6$dose5_selec+c7$dose5_selec+c8$dose5_selec)*100/800

incgen_weib31short$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                               c6$dose1_treat+c7$dose1_treat+c8$dose1_treat)*100/800
incgen_weib31short$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                               c6$dose2_treat+c7$dose2_treat+c8$dose2_treat)*100/800
incgen_weib31short$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                               c6$dose3_treat+c7$dose3_treat+c8$dose3_treat)*100/800
incgen_weib31short$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                               c6$dose4_treat+c7$dose4_treat+c8$dose4_treat)*100/800
incgen_weib31short$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                               c6$dose5_treat+c7$dose5_treat+c8$dose5_treat)*100/800



write.csv(incgen_weib31short,file="incgen_weib31_summaryresults.csv")
