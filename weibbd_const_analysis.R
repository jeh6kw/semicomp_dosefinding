####################Weibull Big Decreasing generation, Constant analysis. Informative Prior and Uninformative Prior


###Inform

for(file in 1:10){
  load(paste("results_weibbdgen_constanalysis_inform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  weibbdinflong <- data.frame(matrix(nrow=5400,ncol=12))
  names(weibbdinflong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  weibbdinflong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  weibbdinflong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  weibbdinflong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      weibbdinflong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      weibbdinflong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      weibbdinflong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      weibbdinflong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      weibbdinflong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      weibbdinflong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      weibbdinflong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      weibbdinflong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      weibbdinflong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  weibbdinflong$chosendoseorder <- as.integer(as.factor(weibbdinflong$chosendose))
  
  
  
  weibbdinfshort <- data.frame(matrix(nrow=54,ncol=16))
  names(weibbdinfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                              "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  weibbdinfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  weibbdinfshort$targetdose <- rep(c(1,3,5),18)
  weibbdinfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- weibbdinflong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(weibbdinflong$duration[((k-1)*100+1):(k*100)],lag=1)
    weibbdinfshort$avgduration[k] <- mean(durationdiff)
    
    weibbdinfshort$dose1_selec[k]<-sum(weibbdinflong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibbdinfshort$dose2_selec[k]<-sum(weibbdinflong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibbdinfshort$dose3_selec[k]<-sum(weibbdinflong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibbdinfshort$dose4_selec[k]<-sum(weibbdinflong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibbdinfshort$dose5_selec[k]<-sum(weibbdinflong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    weibbdinfshort$dose1_treat[k]<- sum(weibbdinflong$numdose1[((k-1)*100+1):(k*100)])*100/sum(weibbdinflong$n[((k-1)*100+1):(k*100)])
    weibbdinfshort$dose2_treat[k]<- sum(weibbdinflong$numdose2[((k-1)*100+1):(k*100)])*100/sum(weibbdinflong$n[((k-1)*100+1):(k*100)])
    weibbdinfshort$dose3_treat[k]<- sum(weibbdinflong$numdose3[((k-1)*100+1):(k*100)])*100/sum(weibbdinflong$n[((k-1)*100+1):(k*100)])
    weibbdinfshort$dose4_treat[k]<- sum(weibbdinflong$numdose4[((k-1)*100+1):(k*100)])*100/sum(weibbdinflong$n[((k-1)*100+1):(k*100)])
    weibbdinfshort$dose5_treat[k]<- sum(weibbdinflong$numdose5[((k-1)*100+1):(k*100)])*100/sum(weibbdinflong$n[((k-1)*100+1):(k*100)])
    
    weibbdinfshort$avgDLT[k] <- mean(weibbdinflong$numDLT[((k-1)*100+1):(k*100)])
    weibbdinfshort$avgDR[k] <- mean(weibbdinflong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(weibbdinfshort,paste("weibbdinf_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("weibbdinf_summaryresults_1.csv")
c2<- read.csv("weibbdinf_summaryresults_2.csv")
c3<- read.csv("weibbdinf_summaryresults_3.csv")
c4<- read.csv("weibbdinf_summaryresults_4.csv")
c5<- read.csv("weibbdinf_summaryresults_5.csv")
c6 <- read.csv("weibbdinf_summaryresults_6.csv")
c7<- read.csv("weibbdinf_summaryresults_7.csv")
c8<- read.csv("weibbdinf_summaryresults_8.csv")
c9<- read.csv("weibbdinf_summaryresults_9.csv")
c10<- read.csv("weibbdinf_summaryresults_10.csv")


weibbdinfshort <- data.frame(matrix(nrow=54,ncol=16))
names(weibbdinfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                            "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
weibbdinfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
weibbdinfshort$targetdose <- rep(c(1,3,5),18)
weibbdinfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

weibbdinfshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                  c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

weibbdinfshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                             c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

weibbdinfshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                            c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

weibbdinfshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                                c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
weibbdinfshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                                c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
weibbdinfshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                                c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
weibbdinfshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                                c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
weibbdinfshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                                c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

weibbdinfshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                 c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
weibbdinfshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                 c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
weibbdinfshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                 c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
weibbdinfshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                 c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
weibbdinfshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                 c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10



write.csv(weibbdinfshort,file="weibbdinf_summaryresults.csv")






###### Uninform

for(file in 1:10){
  load(paste("results_weibbdgen_constanalysis_uninform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  weibbd_uninflong <- data.frame(matrix(nrow=5400,ncol=12))
  names(weibbd_uninflong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  weibbd_uninflong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  weibbd_uninflong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  weibbd_uninflong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      weibbd_uninflong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      weibbd_uninflong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      weibbd_uninflong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      weibbd_uninflong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      weibbd_uninflong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      weibbd_uninflong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      weibbd_uninflong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      weibbd_uninflong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      weibbd_uninflong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  weibbd_uninflong$chosendoseorder <- as.integer(as.factor(weibbd_uninflong$chosendose))
  
  
  
  weibbd_uninfshort <- data.frame(matrix(nrow=54,ncol=16))
  names(weibbd_uninfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                             "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  weibbd_uninfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  weibbd_uninfshort$targetdose <- rep(c(1,3,5),18)
  weibbd_uninfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- weibbd_uninflong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(weibbd_uninflong$duration[((k-1)*100+1):(k*100)],lag=1)
    weibbd_uninfshort$avgduration[k] <- mean(durationdiff)
    
    weibbd_uninfshort$dose1_selec[k]<-sum(weibbd_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibbd_uninfshort$dose2_selec[k]<-sum(weibbd_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibbd_uninfshort$dose3_selec[k]<-sum(weibbd_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibbd_uninfshort$dose4_selec[k]<-sum(weibbd_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibbd_uninfshort$dose5_selec[k]<-sum(weibbd_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    weibbd_uninfshort$dose1_treat[k]<- sum(weibbd_uninflong$numdose1[((k-1)*100+1):(k*100)])*100/sum(weibbd_uninflong$n[((k-1)*100+1):(k*100)])
    weibbd_uninfshort$dose2_treat[k]<- sum(weibbd_uninflong$numdose2[((k-1)*100+1):(k*100)])*100/sum(weibbd_uninflong$n[((k-1)*100+1):(k*100)])
    weibbd_uninfshort$dose3_treat[k]<- sum(weibbd_uninflong$numdose3[((k-1)*100+1):(k*100)])*100/sum(weibbd_uninflong$n[((k-1)*100+1):(k*100)])
    weibbd_uninfshort$dose4_treat[k]<- sum(weibbd_uninflong$numdose4[((k-1)*100+1):(k*100)])*100/sum(weibbd_uninflong$n[((k-1)*100+1):(k*100)])
    weibbd_uninfshort$dose5_treat[k]<- sum(weibbd_uninflong$numdose5[((k-1)*100+1):(k*100)])*100/sum(weibbd_uninflong$n[((k-1)*100+1):(k*100)])
    
    weibbd_uninfshort$avgDLT[k] <- mean(weibbd_uninflong$numDLT[((k-1)*100+1):(k*100)])
    weibbd_uninfshort$avgDR[k] <- mean(weibbd_uninflong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(weibbd_uninfshort,paste("weibbd_uninf_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("weibbd_uninf_summaryresults_1.csv")
c2<- read.csv("weibbd_uninf_summaryresults_2.csv")
c3<- read.csv("weibbd_uninf_summaryresults_3.csv")
c4<- read.csv("weibbd_uninf_summaryresults_4.csv")
c5<- read.csv("weibbd_uninf_summaryresults_5.csv")
c6 <- read.csv("weibbd_uninf_summaryresults_6.csv")
c7<- read.csv("weibbd_uninf_summaryresults_7.csv")
c8<- read.csv("weibbd_uninf_summaryresults_8.csv")
c9<- read.csv("weibbd_uninf_summaryresults_9.csv")
c10<- read.csv("weibbd_uninf_summaryresults_10.csv")


weibbd_uninfshort <- data.frame(matrix(nrow=54,ncol=16))
names(weibbd_uninfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                           "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
weibbd_uninfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
weibbd_uninfshort$targetdose <- rep(c(1,3,5),18)
weibbd_uninfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

weibbd_uninfshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                 c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

weibbd_uninfshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                            c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

weibbd_uninfshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                           c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

weibbd_uninfshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                               c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
weibbd_uninfshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                               c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
weibbd_uninfshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                               c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
weibbd_uninfshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                               c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
weibbd_uninfshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                               c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

weibbd_uninfshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
weibbd_uninfshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
weibbd_uninfshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
weibbd_uninfshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
weibbd_uninfshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10



write.csv(weibbd_uninfshort,file="weibbd_uninf_summaryresults.csv")
