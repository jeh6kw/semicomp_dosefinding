####################All Weibull analysis. C, BD, SD, INC gen. All uninformative priors


###Const gen

for(file in 1:10){
  load(paste("results_constgen_weibanalysis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  const_weiblong <- data.frame(matrix(nrow=5400,ncol=12))
  names(const_weiblong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  const_weiblong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  const_weiblong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  const_weiblong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      const_weiblong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      const_weiblong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      const_weiblong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      const_weiblong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      const_weiblong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      const_weiblong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      const_weiblong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      const_weiblong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      const_weiblong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  const_weiblong$chosendoseorder <- as.integer(as.factor(const_weiblong$chosendose))
  
  
  
  const_weibshort <- data.frame(matrix(nrow=54,ncol=16))
  names(const_weibshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                             "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  const_weibshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  const_weibshort$targetdose <- rep(c(1,3,5),18)
  const_weibshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- const_weiblong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(const_weiblong$duration[((k-1)*100+1):(k*100)],lag=1)
    const_weibshort$avgduration[k] <- mean(durationdiff)
    
    const_weibshort$dose1_selec[k]<-sum(const_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    const_weibshort$dose2_selec[k]<-sum(const_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    const_weibshort$dose3_selec[k]<-sum(const_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    const_weibshort$dose4_selec[k]<-sum(const_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    const_weibshort$dose5_selec[k]<-sum(const_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    const_weibshort$dose1_treat[k]<- sum(const_weiblong$numdose1[((k-1)*100+1):(k*100)])*100/sum(const_weiblong$n[((k-1)*100+1):(k*100)])
    const_weibshort$dose2_treat[k]<- sum(const_weiblong$numdose2[((k-1)*100+1):(k*100)])*100/sum(const_weiblong$n[((k-1)*100+1):(k*100)])
    const_weibshort$dose3_treat[k]<- sum(const_weiblong$numdose3[((k-1)*100+1):(k*100)])*100/sum(const_weiblong$n[((k-1)*100+1):(k*100)])
    const_weibshort$dose4_treat[k]<- sum(const_weiblong$numdose4[((k-1)*100+1):(k*100)])*100/sum(const_weiblong$n[((k-1)*100+1):(k*100)])
    const_weibshort$dose5_treat[k]<- sum(const_weiblong$numdose5[((k-1)*100+1):(k*100)])*100/sum(const_weiblong$n[((k-1)*100+1):(k*100)])
    
    const_weibshort$avgDLT[k] <- mean(const_weiblong$numDLT[((k-1)*100+1):(k*100)])
    const_weibshort$avgDR[k] <- mean(const_weiblong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(const_weibshort,paste("const_weib_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("const_weib_summaryresults_1.csv")
c2<- read.csv("const_weib_summaryresults_2.csv")
c3<- read.csv("const_weib_summaryresults_3.csv")
c4<- read.csv("const_weib_summaryresults_4.csv")
c5<- read.csv("const_weib_summaryresults_5.csv")
c6 <- read.csv("const_weib_summaryresults_6.csv")
c7<- read.csv("const_weib_summaryresults_7.csv")
c8<- read.csv("const_weib_summaryresults_8.csv")
c9<- read.csv("const_weib_summaryresults_9.csv")
c10<- read.csv("const_weib_summaryresults_10.csv")


const_weibshort <- data.frame(matrix(nrow=54,ncol=16))
names(const_weibshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                           "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
const_weibshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
const_weibshort$targetdose <- rep(c(1,3,5),18)
const_weibshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

const_weibshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                 c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

const_weibshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                            c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

const_weibshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                           c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

const_weibshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                               c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
const_weibshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                               c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
const_weibshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                               c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
const_weibshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                               c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
const_weibshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                               c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

const_weibshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
const_weibshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
const_weibshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
const_weibshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
const_weibshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10



write.csv(const_weibshort,file="const_weib_summaryresults.csv")








########################################################################################################
########################################################################################################






### Big Decreasing Gen


for(file in 1:10){
  load(paste("results_weibbdgen_weibanalysis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  bd_weiblong <- data.frame(matrix(nrow=5400,ncol=12))
  names(bd_weiblong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  bd_weiblong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  bd_weiblong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  bd_weiblong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      bd_weiblong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      bd_weiblong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      bd_weiblong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      bd_weiblong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      bd_weiblong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      bd_weiblong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      bd_weiblong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      bd_weiblong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      bd_weiblong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  bd_weiblong$chosendoseorder <- as.integer(as.factor(bd_weiblong$chosendose))
  
  
  
  bd_weibshort <- data.frame(matrix(nrow=54,ncol=16))
  names(bd_weibshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                              "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  bd_weibshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  bd_weibshort$targetdose <- rep(c(1,3,5),18)
  bd_weibshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- bd_weiblong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(bd_weiblong$duration[((k-1)*100+1):(k*100)],lag=1)
    bd_weibshort$avgduration[k] <- mean(durationdiff)
    
    bd_weibshort$dose1_selec[k]<-sum(bd_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    bd_weibshort$dose2_selec[k]<-sum(bd_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    bd_weibshort$dose3_selec[k]<-sum(bd_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    bd_weibshort$dose4_selec[k]<-sum(bd_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    bd_weibshort$dose5_selec[k]<-sum(bd_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    bd_weibshort$dose1_treat[k]<- sum(bd_weiblong$numdose1[((k-1)*100+1):(k*100)])*100/sum(bd_weiblong$n[((k-1)*100+1):(k*100)])
    bd_weibshort$dose2_treat[k]<- sum(bd_weiblong$numdose2[((k-1)*100+1):(k*100)])*100/sum(bd_weiblong$n[((k-1)*100+1):(k*100)])
    bd_weibshort$dose3_treat[k]<- sum(bd_weiblong$numdose3[((k-1)*100+1):(k*100)])*100/sum(bd_weiblong$n[((k-1)*100+1):(k*100)])
    bd_weibshort$dose4_treat[k]<- sum(bd_weiblong$numdose4[((k-1)*100+1):(k*100)])*100/sum(bd_weiblong$n[((k-1)*100+1):(k*100)])
    bd_weibshort$dose5_treat[k]<- sum(bd_weiblong$numdose5[((k-1)*100+1):(k*100)])*100/sum(bd_weiblong$n[((k-1)*100+1):(k*100)])
    
    bd_weibshort$avgDLT[k] <- mean(bd_weiblong$numDLT[((k-1)*100+1):(k*100)])
    bd_weibshort$avgDR[k] <- mean(bd_weiblong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(bd_weibshort,paste("bd_weib_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("bd_weib_summaryresults_1.csv")
c2<- read.csv("bd_weib_summaryresults_2.csv")
c3<- read.csv("bd_weib_summaryresults_3.csv")
c4<- read.csv("bd_weib_summaryresults_4.csv")
c5<- read.csv("bd_weib_summaryresults_5.csv")
c6 <- read.csv("bd_weib_summaryresults_6.csv")
c7<- read.csv("bd_weib_summaryresults_7.csv")
c8<- read.csv("bd_weib_summaryresults_8.csv")
c9<- read.csv("bd_weib_summaryresults_9.csv")
c10<- read.csv("bd_weib_summaryresults_10.csv")


bd_weibshort <- data.frame(matrix(nrow=54,ncol=16))
names(bd_weibshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                            "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
bd_weibshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
bd_weibshort$targetdose <- rep(c(1,3,5),18)
bd_weibshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

bd_weibshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                  c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

bd_weibshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                             c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

bd_weibshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                            c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

bd_weibshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                                c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
bd_weibshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                                c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
bd_weibshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                                c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
bd_weibshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                                c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
bd_weibshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                                c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

bd_weibshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                 c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
bd_weibshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                 c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
bd_weibshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                 c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
bd_weibshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                 c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
bd_weibshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                 c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10



write.csv(bd_weibshort,file="bd_weib_summaryresults.csv")








########################################################################################################
########################################################################################################




### Small Decreasing Gen

for(file in c(1,3:10)){
  load(paste("results_weibsdgen_weibanalysis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  sd_weiblong <- data.frame(matrix(nrow=5400,ncol=12))
  names(sd_weiblong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  sd_weiblong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  sd_weiblong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  sd_weiblong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      sd_weiblong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      sd_weiblong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      sd_weiblong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      sd_weiblong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      sd_weiblong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      sd_weiblong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      sd_weiblong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      sd_weiblong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      sd_weiblong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  sd_weiblong$chosendoseorder <- as.integer(as.factor(sd_weiblong$chosendose))
  
  
  
  sd_weibshort <- data.frame(matrix(nrow=54,ncol=16))
  names(sd_weibshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                           "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  sd_weibshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  sd_weibshort$targetdose <- rep(c(1,3,5),18)
  sd_weibshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- sd_weiblong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(sd_weiblong$duration[((k-1)*100+1):(k*100)],lag=1)
    sd_weibshort$avgduration[k] <- mean(durationdiff)
    
    sd_weibshort$dose1_selec[k]<-sum(sd_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    sd_weibshort$dose2_selec[k]<-sum(sd_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    sd_weibshort$dose3_selec[k]<-sum(sd_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    sd_weibshort$dose4_selec[k]<-sum(sd_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    sd_weibshort$dose5_selec[k]<-sum(sd_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    sd_weibshort$dose1_treat[k]<- sum(sd_weiblong$numdose1[((k-1)*100+1):(k*100)])*100/sum(sd_weiblong$n[((k-1)*100+1):(k*100)])
    sd_weibshort$dose2_treat[k]<- sum(sd_weiblong$numdose2[((k-1)*100+1):(k*100)])*100/sum(sd_weiblong$n[((k-1)*100+1):(k*100)])
    sd_weibshort$dose3_treat[k]<- sum(sd_weiblong$numdose3[((k-1)*100+1):(k*100)])*100/sum(sd_weiblong$n[((k-1)*100+1):(k*100)])
    sd_weibshort$dose4_treat[k]<- sum(sd_weiblong$numdose4[((k-1)*100+1):(k*100)])*100/sum(sd_weiblong$n[((k-1)*100+1):(k*100)])
    sd_weibshort$dose5_treat[k]<- sum(sd_weiblong$numdose5[((k-1)*100+1):(k*100)])*100/sum(sd_weiblong$n[((k-1)*100+1):(k*100)])
    
    sd_weibshort$avgDLT[k] <- mean(sd_weiblong$numDLT[((k-1)*100+1):(k*100)])
    sd_weibshort$avgDR[k] <- mean(sd_weiblong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(sd_weibshort,paste("sd_weib_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("sd_weib_summaryresults_1.csv")
# c2<- read.csv("sd_weib_summaryresults_2.csv")
c3<- read.csv("sd_weib_summaryresults_3.csv")
c4<- read.csv("sd_weib_summaryresults_4.csv")
c5<- read.csv("sd_weib_summaryresults_5.csv")
c6 <- read.csv("sd_weib_summaryresults_6.csv")
c7<- read.csv("sd_weib_summaryresults_7.csv")
c8<- read.csv("sd_weib_summaryresults_8.csv")
c9<- read.csv("sd_weib_summaryresults_9.csv")
c10<- read.csv("sd_weib_summaryresults_10.csv")


sd_weibshort <- data.frame(matrix(nrow=54,ncol=16))
names(sd_weibshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                         "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
sd_weibshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
sd_weibshort$targetdose <- rep(c(1,3,5),18)
sd_weibshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

sd_weibshort$avgduration <- (c1$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                               c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10

sd_weibshort$avgDLT <- (c1$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                          c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)*100/900

sd_weibshort$avgDR <- (c1$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                         c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)*100/900

sd_weibshort$dose1_selec<-(c1$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                             c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)*100/900
sd_weibshort$dose2_selec<-(c1$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                             c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)*100/900
sd_weibshort$dose3_selec<-(c1$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                             c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)*100/900
sd_weibshort$dose4_selec<-(c1$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                             c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)*100/900
sd_weibshort$dose5_selec<-(c1$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                             c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)*100/900

sd_weibshort$dose1_treat<- (c1$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                              c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)*100/900
sd_weibshort$dose2_treat<- (c1$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                              c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)*100/900
sd_weibshort$dose3_treat<- (c1$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                              c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)*100/900
sd_weibshort$dose4_treat<- (c1$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                              c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)*100/900
sd_weibshort$dose5_treat<- (c1$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                              c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)*100/900



write.csv(sd_weibshort,file="sd_weib_summaryresults.csv")







########################################################################################################
########################################################################################################




#### Increasing gen 



for(file in c(1,3:6,8:10)){
  load(paste("results_weibincgen_weibanalysis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  inc_weiblong <- data.frame(matrix(nrow=5400,ncol=12))
  names(inc_weiblong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  inc_weiblong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  inc_weiblong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  inc_weiblong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      inc_weiblong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      inc_weiblong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      inc_weiblong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      inc_weiblong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      inc_weiblong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      inc_weiblong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      inc_weiblong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      inc_weiblong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      inc_weiblong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  inc_weiblong$chosendoseorder <- as.integer(as.factor(inc_weiblong$chosendose))
  
  
  
  inc_weibshort <- data.frame(matrix(nrow=54,ncol=16))
  names(inc_weibshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                           "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  inc_weibshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  inc_weibshort$targetdose <- rep(c(1,3,5),18)
  inc_weibshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- inc_weiblong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(inc_weiblong$duration[((k-1)*100+1):(k*100)],lag=1)
    inc_weibshort$avgduration[k] <- mean(durationdiff)
    
    inc_weibshort$dose1_selec[k]<-sum(inc_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    inc_weibshort$dose2_selec[k]<-sum(inc_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    inc_weibshort$dose3_selec[k]<-sum(inc_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    inc_weibshort$dose4_selec[k]<-sum(inc_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    inc_weibshort$dose5_selec[k]<-sum(inc_weiblong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    inc_weibshort$dose1_treat[k]<- sum(inc_weiblong$numdose1[((k-1)*100+1):(k*100)])*100/sum(inc_weiblong$n[((k-1)*100+1):(k*100)])
    inc_weibshort$dose2_treat[k]<- sum(inc_weiblong$numdose2[((k-1)*100+1):(k*100)])*100/sum(inc_weiblong$n[((k-1)*100+1):(k*100)])
    inc_weibshort$dose3_treat[k]<- sum(inc_weiblong$numdose3[((k-1)*100+1):(k*100)])*100/sum(inc_weiblong$n[((k-1)*100+1):(k*100)])
    inc_weibshort$dose4_treat[k]<- sum(inc_weiblong$numdose4[((k-1)*100+1):(k*100)])*100/sum(inc_weiblong$n[((k-1)*100+1):(k*100)])
    inc_weibshort$dose5_treat[k]<- sum(inc_weiblong$numdose5[((k-1)*100+1):(k*100)])*100/sum(inc_weiblong$n[((k-1)*100+1):(k*100)])
    
    inc_weibshort$avgDLT[k] <- mean(inc_weiblong$numDLT[((k-1)*100+1):(k*100)])
    inc_weibshort$avgDR[k] <- mean(inc_weiblong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(inc_weibshort,paste("inc_weib_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for weibull sd inf
c1 <- read.csv("inc_weib_summaryresults_1.csv")
# c2<- read.csv("inc_weib_summaryresults_2.csv")
c3<- read.csv("inc_weib_summaryresults_3.csv")
c4<- read.csv("inc_weib_summaryresults_4.csv")
c5<- read.csv("inc_weib_summaryresults_5.csv")
c6 <- read.csv("inc_weib_summaryresults_6.csv")
# c7<- read.csv("inc_weib_summaryresults_7.csv")
c8<- read.csv("inc_weib_summaryresults_8.csv")
c9<- read.csv("inc_weib_summaryresults_9.csv")
c10<- read.csv("inc_weib_summaryresults_10.csv")


inc_weibshort <- data.frame(matrix(nrow=54,ncol=16))
names(inc_weibshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                         "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
inc_weibshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
inc_weibshort$targetdose <- rep(c(1,3,5),18)
inc_weibshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

inc_weibshort$avgduration <- (c1$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                               c6$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)*100/800

inc_weibshort$avgDLT <- (c1$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                          c6$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)*100/800

inc_weibshort$avgDR <- (c1$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                         c6$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)*100/800

inc_weibshort$dose1_selec<-(c1$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                             c6$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)*100/800
inc_weibshort$dose2_selec<-(c1$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                             c6$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)*100/800
inc_weibshort$dose3_selec<-(c1$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                             c6$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)*100/800
inc_weibshort$dose4_selec<-(c1$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                             c6$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)*100/800
inc_weibshort$dose5_selec<-(c1$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                             c6$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)*100/800

inc_weibshort$dose1_treat<- (c1$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                              c6$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)*100/800
inc_weibshort$dose2_treat<- (c1$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                              c6$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)*100/800
inc_weibshort$dose3_treat<- (c1$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                              c6$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)*100/800
inc_weibshort$dose4_treat<- (c1$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                              c6$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)*100/800
inc_weibshort$dose5_treat<- (c1$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                              c6$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)*100/800



write.csv(inc_weibshort,file="inc_weib_summaryresults.csv")
