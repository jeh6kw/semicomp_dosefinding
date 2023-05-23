#################### TITECRM MC and CRM MC analysis

###Const gen

for(file in c(1:10)){
  load(paste("results_constgen_titecrm_mcext_",file,".Rdata",sep=''))
  
  cgen_tcrmmclong <- data.frame(matrix(nrow=5400,ncol=13))
  names(cgen_tcrmmclong) <- c("dlttarget","targetdose","chosendose","CRMMCfullinfodose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  cgen_tcrmmclong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  cgen_tcrmmclong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  cgen_tcrmmclong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      cgen_tcrmmclong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      cgen_tcrmmclong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      cgen_tcrmmclong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      cgen_tcrmmclong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      cgen_tcrmmclong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      cgen_tcrmmclong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      cgen_tcrmmclong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      cgen_tcrmmclong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      cgen_tcrmmclong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
      cgen_tcrmmclong$CRMMCfullinfodose[(i-1)*100+j] <- spec$CRMMCfullinfodose[nrow(spec)]
    }
  }
  cgen_tcrmmclong$chosendoseorder <- as.integer(as.factor(cgen_tcrmmclong$chosendose))
  cgen_tcrmmclong$CRMMCfullinfodoseorder <- as.integer(as.factor(cgen_tcrmmclong$CRMMCfullinfodose))
  
  
  cgen_tcrmmcshort <- data.frame(matrix(nrow=54,ncol=21))
  names(cgen_tcrmmcshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                               "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat","dose1_crmmc_selec","dose2_crmmc_selec","dose3_crmmc_selec","dose4_crmmc_selec","dose5_crmmc_selec")
  cgen_tcrmmcshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  cgen_tcrmmcshort$targetdose <- rep(c(1,3,5),18)
  cgen_tcrmmcshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- cgen_tcrmmclong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(cgen_tcrmmclong$duration[((k-1)*100+1):(k*100)],lag=1)
    cgen_tcrmmcshort$avgduration[k] <- mean(durationdiff)
    
    cgen_tcrmmcshort$dose1_selec[k]<-sum(cgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    cgen_tcrmmcshort$dose2_selec[k]<-sum(cgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    cgen_tcrmmcshort$dose3_selec[k]<-sum(cgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    cgen_tcrmmcshort$dose4_selec[k]<-sum(cgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    cgen_tcrmmcshort$dose5_selec[k]<-sum(cgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    cgen_tcrmmcshort$dose1_treat[k]<- sum(cgen_tcrmmclong$numdose1[((k-1)*100+1):(k*100)])*100/sum(cgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    cgen_tcrmmcshort$dose2_treat[k]<- sum(cgen_tcrmmclong$numdose2[((k-1)*100+1):(k*100)])*100/sum(cgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    cgen_tcrmmcshort$dose3_treat[k]<- sum(cgen_tcrmmclong$numdose3[((k-1)*100+1):(k*100)])*100/sum(cgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    cgen_tcrmmcshort$dose4_treat[k]<- sum(cgen_tcrmmclong$numdose4[((k-1)*100+1):(k*100)])*100/sum(cgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    cgen_tcrmmcshort$dose5_treat[k]<- sum(cgen_tcrmmclong$numdose5[((k-1)*100+1):(k*100)])*100/sum(cgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    
    cgen_tcrmmcshort$dose1_crmmc_selec[k]<-sum(cgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    cgen_tcrmmcshort$dose2_crmmc_selec[k]<-sum(cgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    cgen_tcrmmcshort$dose3_crmmc_selec[k]<-sum(cgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    cgen_tcrmmcshort$dose4_crmmc_selec[k]<-sum(cgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    cgen_tcrmmcshort$dose5_crmmc_selec[k]<-sum(cgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    
    cgen_tcrmmcshort$avgDLT[k] <- mean(cgen_tcrmmclong$numDLT[((k-1)*100+1):(k*100)])
    cgen_tcrmmcshort$avgDR[k] <- mean(cgen_tcrmmclong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(cgen_tcrmmcshort,paste("cgen_tcrmmc_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("cgen_tcrmmc_summaryresults_1.csv")
c2<- read.csv("cgen_tcrmmc_summaryresults_2.csv")
c3<- read.csv("cgen_tcrmmc_summaryresults_3.csv")
c4<- read.csv("cgen_tcrmmc_summaryresults_4.csv")
c5<- read.csv("cgen_tcrmmc_summaryresults_5.csv")
c6 <- read.csv("cgen_tcrmmc_summaryresults_6.csv")
c7<- read.csv("cgen_tcrmmc_summaryresults_7.csv")
c8<- read.csv("cgen_tcrmmc_summaryresults_8.csv")
c9<- read.csv("cgen_tcrmmc_summaryresults_9.csv")
c10<- read.csv("cgen_tcrmmc_summaryresults_10.csv")



cgen_tcrmmcshort <- data.frame(matrix(nrow=54,ncol=21))
names(cgen_tcrmmcshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                             "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat","dose1_crmmc_selec","dose2_crmmc_selec","dose3_crmmc_selec","dose4_crmmc_selec","dose5_crmmc_selec")
cgen_tcrmmcshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
cgen_tcrmmcshort$targetdose <- rep(c(1,3,5),18)
cgen_tcrmmcshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

cgen_tcrmmcshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                   c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)*100/1000

cgen_tcrmmcshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                              c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)*100/1000

cgen_tcrmmcshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                             c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)*100/1000

cgen_tcrmmcshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                                 c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)*100/1000
cgen_tcrmmcshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                                 c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)*100/1000
cgen_tcrmmcshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                                 c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)*100/1000
cgen_tcrmmcshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                                 c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)*100/1000
cgen_tcrmmcshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                                 c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)*100/1000

cgen_tcrmmcshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                  c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)*100/1000
cgen_tcrmmcshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                  c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)*100/1000
cgen_tcrmmcshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                  c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)*100/1000
cgen_tcrmmcshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                  c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)*100/1000
cgen_tcrmmcshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                  c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)*100/1000


cgen_tcrmmcshort$dose1_crmmc_selec<-(c1$dose1_crmmc_selec+c2$dose1_crmmc_selec+c3$dose1_crmmc_selec+c4$dose1_crmmc_selec+c5$dose1_crmmc_selec+
                                 c6$dose1_crmmc_selec+c7$dose1_crmmc_selec+c8$dose1_crmmc_selec+c9$dose1_crmmc_selec+c10$dose1_crmmc_selec)*100/1000
cgen_tcrmmcshort$dose2_crmmc_selec<-(c1$dose2_crmmc_selec+c2$dose2_crmmc_selec+c3$dose2_crmmc_selec+c4$dose2_crmmc_selec+c5$dose2_crmmc_selec+
                                 c6$dose2_crmmc_selec+c7$dose2_crmmc_selec+c8$dose2_crmmc_selec+c9$dose2_crmmc_selec+c10$dose2_crmmc_selec)*100/1000
cgen_tcrmmcshort$dose3_crmmc_selec<-(c1$dose3_crmmc_selec+c2$dose3_crmmc_selec+c3$dose3_crmmc_selec+c4$dose3_crmmc_selec+c5$dose3_crmmc_selec+
                                 c6$dose3_crmmc_selec+c7$dose3_crmmc_selec+c8$dose3_crmmc_selec+c9$dose3_crmmc_selec+c10$dose3_crmmc_selec)*100/1000
cgen_tcrmmcshort$dose4_crmmc_selec<-(c1$dose4_crmmc_selec+c2$dose4_crmmc_selec+c3$dose4_crmmc_selec+c4$dose4_crmmc_selec+c5$dose4_crmmc_selec+
                                 c6$dose4_crmmc_selec+c7$dose4_crmmc_selec+c8$dose4_crmmc_selec+c9$dose4_crmmc_selec+c10$dose4_crmmc_selec)*100/1000
cgen_tcrmmcshort$dose5_crmmc_selec<-(c1$dose5_crmmc_selec+c2$dose5_crmmc_selec+c3$dose5_crmmc_selec+c4$dose5_crmmc_selec+c5$dose5_crmmc_selec+
                                 c6$dose5_crmmc_selec+c7$dose5_crmmc_selec+c8$dose5_crmmc_selec+c9$dose5_crmmc_selec+c10$dose5_crmmc_selec)*100/1000

write.csv(cgen_tcrmmcshort,file="cgen_tcrmmc_summaryresults.csv")








########################################################################################################
########################################################################################################






### Big Decreasing Gen

for(file in c(1:10)){
  load(paste("results_weibbdgen_titecrm_mcext_",file,".Rdata",sep=''))
  
  weibbdgen_tcrmmclong <- data.frame(matrix(nrow=5400,ncol=13))
  names(weibbdgen_tcrmmclong) <- c("dlttarget","targetdose","chosendose","CRMMCfullinfodose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  weibbdgen_tcrmmclong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  weibbdgen_tcrmmclong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  weibbdgen_tcrmmclong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      weibbdgen_tcrmmclong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      weibbdgen_tcrmmclong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      weibbdgen_tcrmmclong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      weibbdgen_tcrmmclong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      weibbdgen_tcrmmclong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      weibbdgen_tcrmmclong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      weibbdgen_tcrmmclong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      weibbdgen_tcrmmclong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      weibbdgen_tcrmmclong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
      weibbdgen_tcrmmclong$CRMMCfullinfodose[(i-1)*100+j] <- spec$CRMMCfullinfodose[nrow(spec)]
    }
  }
  weibbdgen_tcrmmclong$chosendoseorder <- as.integer(as.factor(weibbdgen_tcrmmclong$chosendose))
  weibbdgen_tcrmmclong$CRMMCfullinfodoseorder <- as.integer(as.factor(weibbdgen_tcrmmclong$CRMMCfullinfodose))
  
  
  weibbdgen_tcrmmcshort <- data.frame(matrix(nrow=54,ncol=21))
  names(weibbdgen_tcrmmcshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                               "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat","dose1_crmmc_selec","dose2_crmmc_selec","dose3_crmmc_selec","dose4_crmmc_selec","dose5_crmmc_selec")
  weibbdgen_tcrmmcshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  weibbdgen_tcrmmcshort$targetdose <- rep(c(1,3,5),18)
  weibbdgen_tcrmmcshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- weibbdgen_tcrmmclong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(weibbdgen_tcrmmclong$duration[((k-1)*100+1):(k*100)],lag=1)
    weibbdgen_tcrmmcshort$avgduration[k] <- mean(durationdiff)
    
    weibbdgen_tcrmmcshort$dose1_selec[k]<-sum(weibbdgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibbdgen_tcrmmcshort$dose2_selec[k]<-sum(weibbdgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibbdgen_tcrmmcshort$dose3_selec[k]<-sum(weibbdgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibbdgen_tcrmmcshort$dose4_selec[k]<-sum(weibbdgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibbdgen_tcrmmcshort$dose5_selec[k]<-sum(weibbdgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    weibbdgen_tcrmmcshort$dose1_treat[k]<- sum(weibbdgen_tcrmmclong$numdose1[((k-1)*100+1):(k*100)])*100/sum(weibbdgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibbdgen_tcrmmcshort$dose2_treat[k]<- sum(weibbdgen_tcrmmclong$numdose2[((k-1)*100+1):(k*100)])*100/sum(weibbdgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibbdgen_tcrmmcshort$dose3_treat[k]<- sum(weibbdgen_tcrmmclong$numdose3[((k-1)*100+1):(k*100)])*100/sum(weibbdgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibbdgen_tcrmmcshort$dose4_treat[k]<- sum(weibbdgen_tcrmmclong$numdose4[((k-1)*100+1):(k*100)])*100/sum(weibbdgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibbdgen_tcrmmcshort$dose5_treat[k]<- sum(weibbdgen_tcrmmclong$numdose5[((k-1)*100+1):(k*100)])*100/sum(weibbdgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    
    weibbdgen_tcrmmcshort$dose1_crmmc_selec[k]<-sum(weibbdgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibbdgen_tcrmmcshort$dose2_crmmc_selec[k]<-sum(weibbdgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibbdgen_tcrmmcshort$dose3_crmmc_selec[k]<-sum(weibbdgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibbdgen_tcrmmcshort$dose4_crmmc_selec[k]<-sum(weibbdgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibbdgen_tcrmmcshort$dose5_crmmc_selec[k]<-sum(weibbdgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    
    weibbdgen_tcrmmcshort$avgDLT[k] <- mean(weibbdgen_tcrmmclong$numDLT[((k-1)*100+1):(k*100)])
    weibbdgen_tcrmmcshort$avgDR[k] <- mean(weibbdgen_tcrmmclong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(weibbdgen_tcrmmcshort,paste("weibbdgen_tcrmmc_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("weibbdgen_tcrmmc_summaryresults_1.csv")
c2<- read.csv("weibbdgen_tcrmmc_summaryresults_2.csv")
c3<- read.csv("weibbdgen_tcrmmc_summaryresults_3.csv")
c4<- read.csv("weibbdgen_tcrmmc_summaryresults_4.csv")
c5<- read.csv("weibbdgen_tcrmmc_summaryresults_5.csv")
c6 <- read.csv("weibbdgen_tcrmmc_summaryresults_6.csv")
c7<- read.csv("weibbdgen_tcrmmc_summaryresults_7.csv")
c8<- read.csv("weibbdgen_tcrmmc_summaryresults_8.csv")
c9<- read.csv("weibbdgen_tcrmmc_summaryresults_9.csv")
c10<- read.csv("weibbdgen_tcrmmc_summaryresults_10.csv")



weibbdgen_tcrmmcshort <- data.frame(matrix(nrow=54,ncol=21))
names(weibbdgen_tcrmmcshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                             "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat","dose1_crmmc_selec","dose2_crmmc_selec","dose3_crmmc_selec","dose4_crmmc_selec","dose5_crmmc_selec")
weibbdgen_tcrmmcshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
weibbdgen_tcrmmcshort$targetdose <- rep(c(1,3,5),18)
weibbdgen_tcrmmcshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

weibbdgen_tcrmmcshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                   c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)*100/1000

weibbdgen_tcrmmcshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                              c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)*100/1000

weibbdgen_tcrmmcshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                             c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)*100/1000

weibbdgen_tcrmmcshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                                 c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)*100/1000
weibbdgen_tcrmmcshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                                 c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)*100/1000
weibbdgen_tcrmmcshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                                 c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)*100/1000
weibbdgen_tcrmmcshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                                 c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)*100/1000
weibbdgen_tcrmmcshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                                 c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)*100/1000

weibbdgen_tcrmmcshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                  c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)*100/1000
weibbdgen_tcrmmcshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                  c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)*100/1000
weibbdgen_tcrmmcshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                  c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)*100/1000
weibbdgen_tcrmmcshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                  c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)*100/1000
weibbdgen_tcrmmcshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                  c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)*100/1000


weibbdgen_tcrmmcshort$dose1_crmmc_selec<-(c1$dose1_crmmc_selec+c2$dose1_crmmc_selec+c3$dose1_crmmc_selec+c4$dose1_crmmc_selec+c5$dose1_crmmc_selec+
                                       c6$dose1_crmmc_selec+c7$dose1_crmmc_selec+c8$dose1_crmmc_selec+c9$dose1_crmmc_selec+c10$dose1_crmmc_selec)*100/1000
weibbdgen_tcrmmcshort$dose2_crmmc_selec<-(c1$dose2_crmmc_selec+c2$dose2_crmmc_selec+c3$dose2_crmmc_selec+c4$dose2_crmmc_selec+c5$dose2_crmmc_selec+
                                       c6$dose2_crmmc_selec+c7$dose2_crmmc_selec+c8$dose2_crmmc_selec+c9$dose2_crmmc_selec+c10$dose2_crmmc_selec)*100/1000
weibbdgen_tcrmmcshort$dose3_crmmc_selec<-(c1$dose3_crmmc_selec+c2$dose3_crmmc_selec+c3$dose3_crmmc_selec+c4$dose3_crmmc_selec+c5$dose3_crmmc_selec+
                                       c6$dose3_crmmc_selec+c7$dose3_crmmc_selec+c8$dose3_crmmc_selec+c9$dose3_crmmc_selec+c10$dose3_crmmc_selec)*100/1000
weibbdgen_tcrmmcshort$dose4_crmmc_selec<-(c1$dose4_crmmc_selec+c2$dose4_crmmc_selec+c3$dose4_crmmc_selec+c4$dose4_crmmc_selec+c5$dose4_crmmc_selec+
                                       c6$dose4_crmmc_selec+c7$dose4_crmmc_selec+c8$dose4_crmmc_selec+c9$dose4_crmmc_selec+c10$dose4_crmmc_selec)*100/1000
weibbdgen_tcrmmcshort$dose5_crmmc_selec<-(c1$dose5_crmmc_selec+c2$dose5_crmmc_selec+c3$dose5_crmmc_selec+c4$dose5_crmmc_selec+c5$dose5_crmmc_selec+
                                       c6$dose5_crmmc_selec+c7$dose5_crmmc_selec+c8$dose5_crmmc_selec+c9$dose5_crmmc_selec+c10$dose5_crmmc_selec)*100/1000

write.csv(weibbdgen_tcrmmcshort,file="weibbdgen_tcrmmc_summaryresults.csv")






########################################################################################################
########################################################################################################




### Small Decreasing Gen

for(file in c(1:10)){
  load(paste("results_weibsdgen_titecrm_mcext_",file,".Rdata",sep=''))
  
  weibsdgen_tcrmmclong <- data.frame(matrix(nrow=5400,ncol=13))
  names(weibsdgen_tcrmmclong) <- c("dlttarget","targetdose","chosendose","CRMMCfullinfodose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  weibsdgen_tcrmmclong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  weibsdgen_tcrmmclong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  weibsdgen_tcrmmclong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      weibsdgen_tcrmmclong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      weibsdgen_tcrmmclong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      weibsdgen_tcrmmclong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      weibsdgen_tcrmmclong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      weibsdgen_tcrmmclong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      weibsdgen_tcrmmclong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      weibsdgen_tcrmmclong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      weibsdgen_tcrmmclong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      weibsdgen_tcrmmclong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
      weibsdgen_tcrmmclong$CRMMCfullinfodose[(i-1)*100+j] <- spec$CRMMCfullinfodose[nrow(spec)]
    }
  }
  weibsdgen_tcrmmclong$chosendoseorder <- as.integer(as.factor(weibsdgen_tcrmmclong$chosendose))
  weibsdgen_tcrmmclong$CRMMCfullinfodoseorder <- as.integer(as.factor(weibsdgen_tcrmmclong$CRMMCfullinfodose))
  
  
  weibsdgen_tcrmmcshort <- data.frame(matrix(nrow=54,ncol=21))
  names(weibsdgen_tcrmmcshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                               "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat","dose1_crmmc_selec","dose2_crmmc_selec","dose3_crmmc_selec","dose4_crmmc_selec","dose5_crmmc_selec")
  weibsdgen_tcrmmcshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  weibsdgen_tcrmmcshort$targetdose <- rep(c(1,3,5),18)
  weibsdgen_tcrmmcshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- weibsdgen_tcrmmclong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(weibsdgen_tcrmmclong$duration[((k-1)*100+1):(k*100)],lag=1)
    weibsdgen_tcrmmcshort$avgduration[k] <- mean(durationdiff)
    
    weibsdgen_tcrmmcshort$dose1_selec[k]<-sum(weibsdgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibsdgen_tcrmmcshort$dose2_selec[k]<-sum(weibsdgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibsdgen_tcrmmcshort$dose3_selec[k]<-sum(weibsdgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibsdgen_tcrmmcshort$dose4_selec[k]<-sum(weibsdgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibsdgen_tcrmmcshort$dose5_selec[k]<-sum(weibsdgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    weibsdgen_tcrmmcshort$dose1_treat[k]<- sum(weibsdgen_tcrmmclong$numdose1[((k-1)*100+1):(k*100)])*100/sum(weibsdgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibsdgen_tcrmmcshort$dose2_treat[k]<- sum(weibsdgen_tcrmmclong$numdose2[((k-1)*100+1):(k*100)])*100/sum(weibsdgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibsdgen_tcrmmcshort$dose3_treat[k]<- sum(weibsdgen_tcrmmclong$numdose3[((k-1)*100+1):(k*100)])*100/sum(weibsdgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibsdgen_tcrmmcshort$dose4_treat[k]<- sum(weibsdgen_tcrmmclong$numdose4[((k-1)*100+1):(k*100)])*100/sum(weibsdgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibsdgen_tcrmmcshort$dose5_treat[k]<- sum(weibsdgen_tcrmmclong$numdose5[((k-1)*100+1):(k*100)])*100/sum(weibsdgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    
    weibsdgen_tcrmmcshort$dose1_crmmc_selec[k]<-sum(weibsdgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibsdgen_tcrmmcshort$dose2_crmmc_selec[k]<-sum(weibsdgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibsdgen_tcrmmcshort$dose3_crmmc_selec[k]<-sum(weibsdgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibsdgen_tcrmmcshort$dose4_crmmc_selec[k]<-sum(weibsdgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibsdgen_tcrmmcshort$dose5_crmmc_selec[k]<-sum(weibsdgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    
    weibsdgen_tcrmmcshort$avgDLT[k] <- mean(weibsdgen_tcrmmclong$numDLT[((k-1)*100+1):(k*100)])
    weibsdgen_tcrmmcshort$avgDR[k] <- mean(weibsdgen_tcrmmclong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(weibsdgen_tcrmmcshort,paste("weibsdgen_tcrmmc_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("weibsdgen_tcrmmc_summaryresults_1.csv")
c2<- read.csv("weibsdgen_tcrmmc_summaryresults_2.csv")
c3<- read.csv("weibsdgen_tcrmmc_summaryresults_3.csv")
c4<- read.csv("weibsdgen_tcrmmc_summaryresults_4.csv")
c5<- read.csv("weibsdgen_tcrmmc_summaryresults_5.csv")
c6 <- read.csv("weibsdgen_tcrmmc_summaryresults_6.csv")
c7<- read.csv("weibsdgen_tcrmmc_summaryresults_7.csv")
c8<- read.csv("weibsdgen_tcrmmc_summaryresults_8.csv")
c9<- read.csv("weibsdgen_tcrmmc_summaryresults_9.csv")
c10<- read.csv("weibsdgen_tcrmmc_summaryresults_10.csv")



weibsdgen_tcrmmcshort <- data.frame(matrix(nrow=54,ncol=21))
names(weibsdgen_tcrmmcshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                             "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat","dose1_crmmc_selec","dose2_crmmc_selec","dose3_crmmc_selec","dose4_crmmc_selec","dose5_crmmc_selec")
weibsdgen_tcrmmcshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
weibsdgen_tcrmmcshort$targetdose <- rep(c(1,3,5),18)
weibsdgen_tcrmmcshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

weibsdgen_tcrmmcshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                   c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)*100/1000

weibsdgen_tcrmmcshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                              c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)*100/1000

weibsdgen_tcrmmcshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                             c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)*100/1000

weibsdgen_tcrmmcshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                                 c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)*100/1000
weibsdgen_tcrmmcshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                                 c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)*100/1000
weibsdgen_tcrmmcshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                                 c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)*100/1000
weibsdgen_tcrmmcshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                                 c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)*100/1000
weibsdgen_tcrmmcshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                                 c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)*100/1000

weibsdgen_tcrmmcshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                  c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)*100/1000
weibsdgen_tcrmmcshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                  c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)*100/1000
weibsdgen_tcrmmcshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                  c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)*100/1000
weibsdgen_tcrmmcshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                  c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)*100/1000
weibsdgen_tcrmmcshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                  c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)*100/1000


weibsdgen_tcrmmcshort$dose1_crmmc_selec<-(c1$dose1_crmmc_selec+c2$dose1_crmmc_selec+c3$dose1_crmmc_selec+c4$dose1_crmmc_selec+c5$dose1_crmmc_selec+
                                       c6$dose1_crmmc_selec+c7$dose1_crmmc_selec+c8$dose1_crmmc_selec+c9$dose1_crmmc_selec+c10$dose1_crmmc_selec)*100/1000
weibsdgen_tcrmmcshort$dose2_crmmc_selec<-(c1$dose2_crmmc_selec+c2$dose2_crmmc_selec+c3$dose2_crmmc_selec+c4$dose2_crmmc_selec+c5$dose2_crmmc_selec+
                                       c6$dose2_crmmc_selec+c7$dose2_crmmc_selec+c8$dose2_crmmc_selec+c9$dose2_crmmc_selec+c10$dose2_crmmc_selec)*100/1000
weibsdgen_tcrmmcshort$dose3_crmmc_selec<-(c1$dose3_crmmc_selec+c2$dose3_crmmc_selec+c3$dose3_crmmc_selec+c4$dose3_crmmc_selec+c5$dose3_crmmc_selec+
                                       c6$dose3_crmmc_selec+c7$dose3_crmmc_selec+c8$dose3_crmmc_selec+c9$dose3_crmmc_selec+c10$dose3_crmmc_selec)*100/1000
weibsdgen_tcrmmcshort$dose4_crmmc_selec<-(c1$dose4_crmmc_selec+c2$dose4_crmmc_selec+c3$dose4_crmmc_selec+c4$dose4_crmmc_selec+c5$dose4_crmmc_selec+
                                       c6$dose4_crmmc_selec+c7$dose4_crmmc_selec+c8$dose4_crmmc_selec+c9$dose4_crmmc_selec+c10$dose4_crmmc_selec)*100/1000
weibsdgen_tcrmmcshort$dose5_crmmc_selec<-(c1$dose5_crmmc_selec+c2$dose5_crmmc_selec+c3$dose5_crmmc_selec+c4$dose5_crmmc_selec+c5$dose5_crmmc_selec+
                                       c6$dose5_crmmc_selec+c7$dose5_crmmc_selec+c8$dose5_crmmc_selec+c9$dose5_crmmc_selec+c10$dose5_crmmc_selec)*100/1000

write.csv(weibsdgen_tcrmmcshort,file="weibsdgen_tcrmmc_summaryresults.csv")








########################################################################################################
########################################################################################################




#### Increasing gen 

for(file in c(1:10)){
  load(paste("results_weibincgen_titecrm_mcext_",file,".Rdata",sep=''))
  
  weibincgen_tcrmmclong <- data.frame(matrix(nrow=5400,ncol=13))
  names(weibincgen_tcrmmclong) <- c("dlttarget","targetdose","chosendose","CRMMCfullinfodose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  weibincgen_tcrmmclong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  weibincgen_tcrmmclong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  weibincgen_tcrmmclong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      weibincgen_tcrmmclong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      weibincgen_tcrmmclong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      weibincgen_tcrmmclong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      weibincgen_tcrmmclong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      weibincgen_tcrmmclong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      weibincgen_tcrmmclong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      weibincgen_tcrmmclong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      weibincgen_tcrmmclong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      weibincgen_tcrmmclong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
      weibincgen_tcrmmclong$CRMMCfullinfodose[(i-1)*100+j] <- spec$CRMMCfullinfodose[nrow(spec)]
    }
  }
  weibincgen_tcrmmclong$chosendoseorder <- as.integer(as.factor(weibincgen_tcrmmclong$chosendose))
  weibincgen_tcrmmclong$CRMMCfullinfodoseorder <- as.integer(as.factor(weibincgen_tcrmmclong$CRMMCfullinfodose))
  
  
  weibincgen_tcrmmcshort <- data.frame(matrix(nrow=54,ncol=21))
  names(weibincgen_tcrmmcshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                               "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat","dose1_crmmc_selec","dose2_crmmc_selec","dose3_crmmc_selec","dose4_crmmc_selec","dose5_crmmc_selec")
  weibincgen_tcrmmcshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  weibincgen_tcrmmcshort$targetdose <- rep(c(1,3,5),18)
  weibincgen_tcrmmcshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- weibincgen_tcrmmclong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(weibincgen_tcrmmclong$duration[((k-1)*100+1):(k*100)],lag=1)
    weibincgen_tcrmmcshort$avgduration[k] <- mean(durationdiff)
    
    weibincgen_tcrmmcshort$dose1_selec[k]<-sum(weibincgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibincgen_tcrmmcshort$dose2_selec[k]<-sum(weibincgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibincgen_tcrmmcshort$dose3_selec[k]<-sum(weibincgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibincgen_tcrmmcshort$dose4_selec[k]<-sum(weibincgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibincgen_tcrmmcshort$dose5_selec[k]<-sum(weibincgen_tcrmmclong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    weibincgen_tcrmmcshort$dose1_treat[k]<- sum(weibincgen_tcrmmclong$numdose1[((k-1)*100+1):(k*100)])*100/sum(weibincgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibincgen_tcrmmcshort$dose2_treat[k]<- sum(weibincgen_tcrmmclong$numdose2[((k-1)*100+1):(k*100)])*100/sum(weibincgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibincgen_tcrmmcshort$dose3_treat[k]<- sum(weibincgen_tcrmmclong$numdose3[((k-1)*100+1):(k*100)])*100/sum(weibincgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibincgen_tcrmmcshort$dose4_treat[k]<- sum(weibincgen_tcrmmclong$numdose4[((k-1)*100+1):(k*100)])*100/sum(weibincgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    weibincgen_tcrmmcshort$dose5_treat[k]<- sum(weibincgen_tcrmmclong$numdose5[((k-1)*100+1):(k*100)])*100/sum(weibincgen_tcrmmclong$n[((k-1)*100+1):(k*100)])
    
    weibincgen_tcrmmcshort$dose1_crmmc_selec[k]<-sum(weibincgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    weibincgen_tcrmmcshort$dose2_crmmc_selec[k]<-sum(weibincgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    weibincgen_tcrmmcshort$dose3_crmmc_selec[k]<-sum(weibincgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    weibincgen_tcrmmcshort$dose4_crmmc_selec[k]<-sum(weibincgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    weibincgen_tcrmmcshort$dose5_crmmc_selec[k]<-sum(weibincgen_tcrmmclong$CRMMCfullinfodoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    
    weibincgen_tcrmmcshort$avgDLT[k] <- mean(weibincgen_tcrmmclong$numDLT[((k-1)*100+1):(k*100)])
    weibincgen_tcrmmcshort$avgDR[k] <- mean(weibincgen_tcrmmclong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(weibincgen_tcrmmcshort,paste("weibincgen_tcrmmc_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("weibincgen_tcrmmc_summaryresults_1.csv")
c2<- read.csv("weibincgen_tcrmmc_summaryresults_2.csv")
c3<- read.csv("weibincgen_tcrmmc_summaryresults_3.csv")
c4<- read.csv("weibincgen_tcrmmc_summaryresults_4.csv")
c5<- read.csv("weibincgen_tcrmmc_summaryresults_5.csv")
c6 <- read.csv("weibincgen_tcrmmc_summaryresults_6.csv")
c7<- read.csv("weibincgen_tcrmmc_summaryresults_7.csv")
c8<- read.csv("weibincgen_tcrmmc_summaryresults_8.csv")
c9<- read.csv("weibincgen_tcrmmc_summaryresults_9.csv")
c10<- read.csv("weibincgen_tcrmmc_summaryresults_10.csv")



weibincgen_tcrmmcshort <- data.frame(matrix(nrow=54,ncol=21))
names(weibincgen_tcrmmcshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                             "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat","dose1_crmmc_selec","dose2_crmmc_selec","dose3_crmmc_selec","dose4_crmmc_selec","dose5_crmmc_selec")
weibincgen_tcrmmcshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
weibincgen_tcrmmcshort$targetdose <- rep(c(1,3,5),18)
weibincgen_tcrmmcshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

weibincgen_tcrmmcshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                   c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)*100/1000

weibincgen_tcrmmcshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                              c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)*100/1000

weibincgen_tcrmmcshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                             c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)*100/1000

weibincgen_tcrmmcshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                                 c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)*100/1000
weibincgen_tcrmmcshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                                 c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)*100/1000
weibincgen_tcrmmcshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                                 c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)*100/1000
weibincgen_tcrmmcshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                                 c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)*100/1000
weibincgen_tcrmmcshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                                 c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)*100/1000

weibincgen_tcrmmcshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                  c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)*100/1000
weibincgen_tcrmmcshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                  c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)*100/1000
weibincgen_tcrmmcshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                  c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)*100/1000
weibincgen_tcrmmcshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                  c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)*100/1000
weibincgen_tcrmmcshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                  c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)*100/1000


weibincgen_tcrmmcshort$dose1_crmmc_selec<-(c1$dose1_crmmc_selec+c2$dose1_crmmc_selec+c3$dose1_crmmc_selec+c4$dose1_crmmc_selec+c5$dose1_crmmc_selec+
                                       c6$dose1_crmmc_selec+c7$dose1_crmmc_selec+c8$dose1_crmmc_selec+c9$dose1_crmmc_selec+c10$dose1_crmmc_selec)*100/1000
weibincgen_tcrmmcshort$dose2_crmmc_selec<-(c1$dose2_crmmc_selec+c2$dose2_crmmc_selec+c3$dose2_crmmc_selec+c4$dose2_crmmc_selec+c5$dose2_crmmc_selec+
                                       c6$dose2_crmmc_selec+c7$dose2_crmmc_selec+c8$dose2_crmmc_selec+c9$dose2_crmmc_selec+c10$dose2_crmmc_selec)*100/1000
weibincgen_tcrmmcshort$dose3_crmmc_selec<-(c1$dose3_crmmc_selec+c2$dose3_crmmc_selec+c3$dose3_crmmc_selec+c4$dose3_crmmc_selec+c5$dose3_crmmc_selec+
                                       c6$dose3_crmmc_selec+c7$dose3_crmmc_selec+c8$dose3_crmmc_selec+c9$dose3_crmmc_selec+c10$dose3_crmmc_selec)*100/1000
weibincgen_tcrmmcshort$dose4_crmmc_selec<-(c1$dose4_crmmc_selec+c2$dose4_crmmc_selec+c3$dose4_crmmc_selec+c4$dose4_crmmc_selec+c5$dose4_crmmc_selec+
                                       c6$dose4_crmmc_selec+c7$dose4_crmmc_selec+c8$dose4_crmmc_selec+c9$dose4_crmmc_selec+c10$dose4_crmmc_selec)*100/1000
weibincgen_tcrmmcshort$dose5_crmmc_selec<-(c1$dose5_crmmc_selec+c2$dose5_crmmc_selec+c3$dose5_crmmc_selec+c4$dose5_crmmc_selec+c5$dose5_crmmc_selec+
                                       c6$dose5_crmmc_selec+c7$dose5_crmmc_selec+c8$dose5_crmmc_selec+c9$dose5_crmmc_selec+c10$dose5_crmmc_selec)*100/1000

write.csv(weibincgen_tcrmmcshort,file="weibincgen_tcrmmc_summaryresults.csv")





