####################Weibull Small Decreasing generation, Constant analysis. Informative Prior and Uninformative Prior


###Inform

for(file in 1:10){
  load(paste("results_weibsdgen_constanalysis_inform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  weibsdinflong <- data.frame(matrix(nrow=5400,ncol=12))
  names(weibsdinflong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  weibsdinflong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  weibsdinflong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  weibsdinflong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
   
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      weibsdinflong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      weibsdinflong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      weibsdinflong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      weibsdinflong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      weibsdinflong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      weibsdinflong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      weibsdinflong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      weibsdinflong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      weibsdinflong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  weibsdinflong$chosendoseorder <- as.integer(as.factor(weibsdinflong$chosendose))
  
  
  
  weibsdinfshort <- data.frame(matrix(nrow=54,ncol=16))
  names(weibsdinfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                            "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  weibsdinfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  weibsdinfshort$targetdose <- rep(c(1,3,5),18)
  weibsdinfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- weibsdinflong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(weibsdinflong$duration[((k-1)*100+1):(k*100)],lag=1)
    weibsdinfshort$avgduration[k] <- mean(durationdiff)
    
    weibsdinfshort$dose1_selec[k]<-sum(weibsdinflong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibsdinfshort$dose2_selec[k]<-sum(weibsdinflong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibsdinfshort$dose3_selec[k]<-sum(weibsdinflong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibsdinfshort$dose4_selec[k]<-sum(weibsdinflong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibsdinfshort$dose5_selec[k]<-sum(weibsdinflong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    weibsdinfshort$dose1_treat[k]<- sum(weibsdinflong$numdose1[((k-1)*100+1):(k*100)])*100/sum(weibsdinflong$n[((k-1)*100+1):(k*100)])
    weibsdinfshort$dose2_treat[k]<- sum(weibsdinflong$numdose2[((k-1)*100+1):(k*100)])*100/sum(weibsdinflong$n[((k-1)*100+1):(k*100)])
    weibsdinfshort$dose3_treat[k]<- sum(weibsdinflong$numdose3[((k-1)*100+1):(k*100)])*100/sum(weibsdinflong$n[((k-1)*100+1):(k*100)])
    weibsdinfshort$dose4_treat[k]<- sum(weibsdinflong$numdose4[((k-1)*100+1):(k*100)])*100/sum(weibsdinflong$n[((k-1)*100+1):(k*100)])
    weibsdinfshort$dose5_treat[k]<- sum(weibsdinflong$numdose5[((k-1)*100+1):(k*100)])*100/sum(weibsdinflong$n[((k-1)*100+1):(k*100)])
    
    weibsdinfshort$avgDLT[k] <- mean(weibsdinflong$numDLT[((k-1)*100+1):(k*100)])
    weibsdinfshort$avgDR[k] <- mean(weibsdinflong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(weibsdinfshort,paste("weibsdinf_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("weibsdinf_summaryresults_1.csv")
c2<- read.csv("weibsdinf_summaryresults_2.csv")
c3<- read.csv("weibsdinf_summaryresults_3.csv")
c4<- read.csv("weibsdinf_summaryresults_4.csv")
c5<- read.csv("weibsdinf_summaryresults_5.csv")
c6 <- read.csv("weibsdinf_summaryresults_6.csv")
c7<- read.csv("weibsdinf_summaryresults_7.csv")
c8<- read.csv("weibsdinf_summaryresults_8.csv")
c9<- read.csv("weibsdinf_summaryresults_9.csv")
c10<- read.csv("weibsdinf_summaryresults_10.csv")


weibsdinfshort <- data.frame(matrix(nrow=54,ncol=16))
names(weibsdinfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                          "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
weibsdinfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
weibsdinfshort$targetdose <- rep(c(1,3,5),18)
weibsdinfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

weibsdinfshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

weibsdinfshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                                 c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

weibsdinfshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                            c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

weibsdinfshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                              c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
weibsdinfshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                              c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
weibsdinfshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                              c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
weibsdinfshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                              c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
weibsdinfshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                              c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

weibsdinfshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                               c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
weibsdinfshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                               c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
weibsdinfshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                               c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
weibsdinfshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                               c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
weibsdinfshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                               c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10



write.csv(weibsdinfshort,file="weibsdinf_summaryresults.csv")










###Uninform Normal

####################Weibull Small Decreasing generation, Constant analysis. Informative Prior and Uninformative Prior


for(file in 1:10){
  load(paste("results_weibsdgen_constanalysis_uninform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  weibsd_uninflong <- data.frame(matrix(nrow=5400,ncol=12))
  names(weibsd_uninflong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  weibsd_uninflong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  weibsd_uninflong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  weibsd_uninflong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      weibsd_uninflong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      weibsd_uninflong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      weibsd_uninflong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      weibsd_uninflong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      weibsd_uninflong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      weibsd_uninflong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      weibsd_uninflong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      weibsd_uninflong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      weibsd_uninflong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  weibsd_uninflong$chosendoseorder <- as.integer(as.factor(weibsd_uninflong$chosendose))
  
  
  
  weibsd_uninfshort <- data.frame(matrix(nrow=54,ncol=16))
  names(weibsd_uninfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                             "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  weibsd_uninfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  weibsd_uninfshort$targetdose <- rep(c(1,3,5),18)
  weibsd_uninfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- weibsd_uninflong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(weibsd_uninflong$duration[((k-1)*100+1):(k*100)],lag=1)
    weibsd_uninfshort$avgduration[k] <- mean(durationdiff)
    
    weibsd_uninfshort$dose1_selec[k]<-sum(weibsd_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibsd_uninfshort$dose2_selec[k]<-sum(weibsd_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibsd_uninfshort$dose3_selec[k]<-sum(weibsd_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibsd_uninfshort$dose4_selec[k]<-sum(weibsd_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibsd_uninfshort$dose5_selec[k]<-sum(weibsd_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    weibsd_uninfshort$dose1_treat[k]<- sum(weibsd_uninflong$numdose1[((k-1)*100+1):(k*100)])*100/sum(weibsd_uninflong$n[((k-1)*100+1):(k*100)])
    weibsd_uninfshort$dose2_treat[k]<- sum(weibsd_uninflong$numdose2[((k-1)*100+1):(k*100)])*100/sum(weibsd_uninflong$n[((k-1)*100+1):(k*100)])
    weibsd_uninfshort$dose3_treat[k]<- sum(weibsd_uninflong$numdose3[((k-1)*100+1):(k*100)])*100/sum(weibsd_uninflong$n[((k-1)*100+1):(k*100)])
    weibsd_uninfshort$dose4_treat[k]<- sum(weibsd_uninflong$numdose4[((k-1)*100+1):(k*100)])*100/sum(weibsd_uninflong$n[((k-1)*100+1):(k*100)])
    weibsd_uninfshort$dose5_treat[k]<- sum(weibsd_uninflong$numdose5[((k-1)*100+1):(k*100)])*100/sum(weibsd_uninflong$n[((k-1)*100+1):(k*100)])
    
    weibsd_uninfshort$avgDLT[k] <- mean(weibsd_uninflong$numDLT[((k-1)*100+1):(k*100)])
    weibsd_uninfshort$avgDR[k] <- mean(weibsd_uninflong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(weibsd_uninfshort,paste("weibsd_uninf_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("weibsd_uninf_summaryresults_1.csv")
c2<- read.csv("weibsd_uninf_summaryresults_2.csv")
c3<- read.csv("weibsd_uninf_summaryresults_3.csv")
c4<- read.csv("weibsd_uninf_summaryresults_4.csv")
c5<- read.csv("weibsd_uninf_summaryresults_5.csv")
c6 <- read.csv("weibsd_uninf_summaryresults_6.csv")
c7<- read.csv("weibsd_uninf_summaryresults_7.csv")
c8<- read.csv("weibsd_uninf_summaryresults_8.csv")
c9<- read.csv("weibsd_uninf_summaryresults_9.csv")
c10<- read.csv("weibsd_uninf_summaryresults_10.csv")


weibsd_uninfshort <- data.frame(matrix(nrow=54,ncol=16))
names(weibsd_uninfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                           "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
weibsd_uninfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
weibsd_uninfshort$targetdose <- rep(c(1,3,5),18)
weibsd_uninfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

weibsd_uninfshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                 c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

weibsd_uninfshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                            c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

weibsd_uninfshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                           c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

weibsd_uninfshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                               c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
weibsd_uninfshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                               c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
weibsd_uninfshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                               c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
weibsd_uninfshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                               c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
weibsd_uninfshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                               c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

weibsd_uninfshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
weibsd_uninfshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
weibsd_uninfshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
weibsd_uninfshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
weibsd_uninfshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10



write.csv(weibsd_uninfshort,file="weibsd_uninf_summaryresults.csv")






################ Uniform prior

for(file in 1:10){
  load(paste("results_weibsdgen_constanalysis_uniformprior_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  weibsd_uniformpriorlong <- data.frame(matrix(nrow=5400,ncol=12))
  names(weibsd_uniformpriorlong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  weibsd_uniformpriorlong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  weibsd_uniformpriorlong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  weibsd_uniformpriorlong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      weibsd_uniformpriorlong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      weibsd_uniformpriorlong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      weibsd_uniformpriorlong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      weibsd_uniformpriorlong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      weibsd_uniformpriorlong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      weibsd_uniformpriorlong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      weibsd_uniformpriorlong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      weibsd_uniformpriorlong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      weibsd_uniformpriorlong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  weibsd_uniformpriorlong$chosendoseorder <- as.integer(as.factor(weibsd_uniformpriorlong$chosendose))
  
  
  
  weibsd_uniformpriorshort <- data.frame(matrix(nrow=54,ncol=16))
  names(weibsd_uniformpriorshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                                "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  weibsd_uniformpriorshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  weibsd_uniformpriorshort$targetdose <- rep(c(1,3,5),18)
  weibsd_uniformpriorshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- weibsd_uniformpriorlong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(weibsd_uniformpriorlong$duration[((k-1)*100+1):(k*100)],lag=1)
    weibsd_uniformpriorshort$avgduration[k] <- mean(durationdiff)
    
    weibsd_uniformpriorshort$dose1_selec[k]<-sum(weibsd_uniformpriorlong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibsd_uniformpriorshort$dose2_selec[k]<-sum(weibsd_uniformpriorlong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibsd_uniformpriorshort$dose3_selec[k]<-sum(weibsd_uniformpriorlong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibsd_uniformpriorshort$dose4_selec[k]<-sum(weibsd_uniformpriorlong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibsd_uniformpriorshort$dose5_selec[k]<-sum(weibsd_uniformpriorlong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    weibsd_uniformpriorshort$dose1_treat[k]<- sum(weibsd_uniformpriorlong$numdose1[((k-1)*100+1):(k*100)])*100/sum(weibsd_uniformpriorlong$n[((k-1)*100+1):(k*100)])
    weibsd_uniformpriorshort$dose2_treat[k]<- sum(weibsd_uniformpriorlong$numdose2[((k-1)*100+1):(k*100)])*100/sum(weibsd_uniformpriorlong$n[((k-1)*100+1):(k*100)])
    weibsd_uniformpriorshort$dose3_treat[k]<- sum(weibsd_uniformpriorlong$numdose3[((k-1)*100+1):(k*100)])*100/sum(weibsd_uniformpriorlong$n[((k-1)*100+1):(k*100)])
    weibsd_uniformpriorshort$dose4_treat[k]<- sum(weibsd_uniformpriorlong$numdose4[((k-1)*100+1):(k*100)])*100/sum(weibsd_uniformpriorlong$n[((k-1)*100+1):(k*100)])
    weibsd_uniformpriorshort$dose5_treat[k]<- sum(weibsd_uniformpriorlong$numdose5[((k-1)*100+1):(k*100)])*100/sum(weibsd_uniformpriorlong$n[((k-1)*100+1):(k*100)])
    
    weibsd_uniformpriorshort$avgDLT[k] <- mean(weibsd_uniformpriorlong$numDLT[((k-1)*100+1):(k*100)])
    weibsd_uniformpriorshort$avgDR[k] <- mean(weibsd_uniformpriorlong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(weibsd_uniformpriorshort,paste("weibsd_uniformprior_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("weibsd_uniformprior_summaryresults_1.csv")
c2<- read.csv("weibsd_uniformprior_summaryresults_2.csv")
c3<- read.csv("weibsd_uniformprior_summaryresults_3.csv")
c4<- read.csv("weibsd_uniformprior_summaryresults_4.csv")
c5<- read.csv("weibsd_uniformprior_summaryresults_5.csv")
c6 <- read.csv("weibsd_uniformprior_summaryresults_6.csv")
c7<- read.csv("weibsd_uniformprior_summaryresults_7.csv")
c8<- read.csv("weibsd_uniformprior_summaryresults_8.csv")
c9<- read.csv("weibsd_uniformprior_summaryresults_9.csv")
c10<- read.csv("weibsd_uniformprior_summaryresults_10.csv")


weibsd_uniformpriorshort <- data.frame(matrix(nrow=54,ncol=16))
names(weibsd_uniformpriorshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                              "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
weibsd_uniformpriorshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
weibsd_uniformpriorshort$targetdose <- rep(c(1,3,5),18)
weibsd_uniformpriorshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

weibsd_uniformpriorshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                    c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

weibsd_uniformpriorshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                               c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

weibsd_uniformpriorshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                              c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

weibsd_uniformpriorshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                                  c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
weibsd_uniformpriorshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                                  c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
weibsd_uniformpriorshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                                  c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
weibsd_uniformpriorshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                                  c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
weibsd_uniformpriorshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                                  c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

weibsd_uniformpriorshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                   c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
weibsd_uniformpriorshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                   c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
weibsd_uniformpriorshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                   c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
weibsd_uniformpriorshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                   c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
weibsd_uniformpriorshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                   c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10



write.csv(weibsd_uniformpriorshort,file="weibsd_uniformprior_summaryresults.csv")
