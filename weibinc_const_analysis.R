####################Weibull Small Decreasing generation, Constant analysis. Informative Prior and Uninformative Prior


###Inform

for(file in 1:10){
  load(paste("results_weibincgen_constanalysis_inform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  weibincinflong <- data.frame(matrix(nrow=5400,ncol=12))
  names(weibincinflong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  weibincinflong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  weibincinflong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  weibincinflong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      weibincinflong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      weibincinflong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      weibincinflong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      weibincinflong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      weibincinflong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      weibincinflong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      weibincinflong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      weibincinflong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      weibincinflong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  weibincinflong$chosendoseorder <- as.integer(as.factor(weibincinflong$chosendose))
  
  
  
  weibincinfshort <- data.frame(matrix(nrow=54,ncol=16))
  names(weibincinfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                             "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  weibincinfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  weibincinfshort$targetdose <- rep(c(1,3,5),18)
  weibincinfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- weibincinflong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(weibincinflong$duration[((k-1)*100+1):(k*100)],lag=1)
    weibincinfshort$avgduration[k] <- mean(durationdiff)
    
    weibincinfshort$dose1_selec[k]<-sum(weibincinflong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibincinfshort$dose2_selec[k]<-sum(weibincinflong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibincinfshort$dose3_selec[k]<-sum(weibincinflong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibincinfshort$dose4_selec[k]<-sum(weibincinflong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibincinfshort$dose5_selec[k]<-sum(weibincinflong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    weibincinfshort$dose1_treat[k]<- sum(weibincinflong$numdose1[((k-1)*100+1):(k*100)])*100/sum(weibincinflong$n[((k-1)*100+1):(k*100)])
    weibincinfshort$dose2_treat[k]<- sum(weibincinflong$numdose2[((k-1)*100+1):(k*100)])*100/sum(weibincinflong$n[((k-1)*100+1):(k*100)])
    weibincinfshort$dose3_treat[k]<- sum(weibincinflong$numdose3[((k-1)*100+1):(k*100)])*100/sum(weibincinflong$n[((k-1)*100+1):(k*100)])
    weibincinfshort$dose4_treat[k]<- sum(weibincinflong$numdose4[((k-1)*100+1):(k*100)])*100/sum(weibincinflong$n[((k-1)*100+1):(k*100)])
    weibincinfshort$dose5_treat[k]<- sum(weibincinflong$numdose5[((k-1)*100+1):(k*100)])*100/sum(weibincinflong$n[((k-1)*100+1):(k*100)])
    
    weibincinfshort$avgDLT[k] <- mean(weibincinflong$numDLT[((k-1)*100+1):(k*100)])
    weibincinfshort$avgDR[k] <- mean(weibincinflong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(weibincinfshort,paste("weibincinf_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("weibincinf_summaryresults_1.csv")
c2<- read.csv("weibincinf_summaryresults_2.csv")
c3<- read.csv("weibincinf_summaryresults_3.csv")
c4<- read.csv("weibincinf_summaryresults_4.csv")
c5<- read.csv("weibincinf_summaryresults_5.csv")
c6 <- read.csv("weibincinf_summaryresults_6.csv")
c7<- read.csv("weibincinf_summaryresults_7.csv")
c8<- read.csv("weibincinf_summaryresults_8.csv")
c9<- read.csv("weibincinf_summaryresults_9.csv")
c10<- read.csv("weibincinf_summaryresults_10.csv")


weibincinfshort <- data.frame(matrix(nrow=54,ncol=16))
names(weibincinfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                           "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
weibincinfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
weibincinfshort$targetdose <- rep(c(1,3,5),18)
weibincinfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

weibincinfshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                 c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

weibincinfshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                            c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

weibincinfshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                           c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

weibincinfshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                               c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
weibincinfshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                               c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
weibincinfshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                               c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
weibincinfshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                               c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
weibincinfshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                               c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

weibincinfshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
weibincinfshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
weibincinfshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
weibincinfshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
weibincinfshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10



write.csv(weibincinfshort,file="weibincinf_summaryresults.csv")






###Uninform

for(file in 1:10){
  load(paste("results_weibincgen_constanalysis_uninform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  weibinc_uninflong <- data.frame(matrix(nrow=5400,ncol=12))
  names(weibinc_uninflong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  weibinc_uninflong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  weibinc_uninflong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  weibinc_uninflong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      weibinc_uninflong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      weibinc_uninflong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      weibinc_uninflong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      weibinc_uninflong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      weibinc_uninflong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      weibinc_uninflong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      weibinc_uninflong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      weibinc_uninflong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      weibinc_uninflong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  weibinc_uninflong$chosendoseorder <- as.integer(as.factor(weibinc_uninflong$chosendose))
  
  
  
  weibinc_uninfshort <- data.frame(matrix(nrow=54,ncol=16))
  names(weibinc_uninfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                              "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  weibinc_uninfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  weibinc_uninfshort$targetdose <- rep(c(1,3,5),18)
  weibinc_uninfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- weibinc_uninflong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(weibinc_uninflong$duration[((k-1)*100+1):(k*100)],lag=1)
    weibinc_uninfshort$avgduration[k] <- mean(durationdiff)
    
    weibinc_uninfshort$dose1_selec[k]<-sum(weibinc_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibinc_uninfshort$dose2_selec[k]<-sum(weibinc_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibinc_uninfshort$dose3_selec[k]<-sum(weibinc_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibinc_uninfshort$dose4_selec[k]<-sum(weibinc_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibinc_uninfshort$dose5_selec[k]<-sum(weibinc_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    weibinc_uninfshort$dose1_treat[k]<- sum(weibinc_uninflong$numdose1[((k-1)*100+1):(k*100)])*100/sum(weibinc_uninflong$n[((k-1)*100+1):(k*100)])
    weibinc_uninfshort$dose2_treat[k]<- sum(weibinc_uninflong$numdose2[((k-1)*100+1):(k*100)])*100/sum(weibinc_uninflong$n[((k-1)*100+1):(k*100)])
    weibinc_uninfshort$dose3_treat[k]<- sum(weibinc_uninflong$numdose3[((k-1)*100+1):(k*100)])*100/sum(weibinc_uninflong$n[((k-1)*100+1):(k*100)])
    weibinc_uninfshort$dose4_treat[k]<- sum(weibinc_uninflong$numdose4[((k-1)*100+1):(k*100)])*100/sum(weibinc_uninflong$n[((k-1)*100+1):(k*100)])
    weibinc_uninfshort$dose5_treat[k]<- sum(weibinc_uninflong$numdose5[((k-1)*100+1):(k*100)])*100/sum(weibinc_uninflong$n[((k-1)*100+1):(k*100)])
    
    weibinc_uninfshort$avgDLT[k] <- mean(weibinc_uninflong$numDLT[((k-1)*100+1):(k*100)])
    weibinc_uninfshort$avgDR[k] <- mean(weibinc_uninflong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(weibinc_uninfshort,paste("weibinc_uninf_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("weibinc_uninf_summaryresults_1.csv")
c2<- read.csv("weibinc_uninf_summaryresults_2.csv")
c3<- read.csv("weibinc_uninf_summaryresults_3.csv")
c4<- read.csv("weibinc_uninf_summaryresults_4.csv")
c5<- read.csv("weibinc_uninf_summaryresults_5.csv")
c6 <- read.csv("weibinc_uninf_summaryresults_6.csv")
c7<- read.csv("weibinc_uninf_summaryresults_7.csv")
c8<- read.csv("weibinc_uninf_summaryresults_8.csv")
c9<- read.csv("weibinc_uninf_summaryresults_9.csv")
c10<- read.csv("weibinc_uninf_summaryresults_10.csv")


weibinc_uninfshort <- data.frame(matrix(nrow=54,ncol=16))
names(weibinc_uninfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                            "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
weibinc_uninfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
weibinc_uninfshort$targetdose <- rep(c(1,3,5),18)
weibinc_uninfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

weibinc_uninfshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                  c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

weibinc_uninfshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                             c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

weibinc_uninfshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                            c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

weibinc_uninfshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                                c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
weibinc_uninfshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                                c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
weibinc_uninfshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                                c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
weibinc_uninfshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                                c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
weibinc_uninfshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                                c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

weibinc_uninfshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                 c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
weibinc_uninfshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                 c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
weibinc_uninfshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                 c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
weibinc_uninfshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                 c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
weibinc_uninfshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                 c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10



write.csv(weibinc_uninfshort,file="weibinc_uninf_summaryresults.csv")






