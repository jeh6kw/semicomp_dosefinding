####################Constant generation, Constant analysis. 


###Inform Normal prior

for(file in 1:10){
  load(paste("results_constgen_constanalysis_inform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  constinflong <- data.frame(matrix(nrow=5400,ncol=12))
  names(constinflong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  constinflong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  constinflong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  constinflong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      constinflong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      constinflong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      constinflong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      constinflong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      constinflong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      constinflong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      constinflong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      constinflong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      constinflong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  constinflong$chosendoseorder <- as.integer(as.factor(constinflong$chosendose))
  
  
  
  constinfshort <- data.frame(matrix(nrow=54,ncol=16))
  names(constinfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                            "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  constinfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  constinfshort$targetdose <- rep(c(1,3,5),18)
  constinfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- constinflong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(constinflong$duration[((k-1)*100+1):(k*100)],lag=1)
    constinfshort$avgduration[k] <- mean(durationdiff)
    
    constinfshort$dose1_selec[k]<-sum(constinflong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    constinfshort$dose2_selec[k]<-sum(constinflong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    constinfshort$dose3_selec[k]<-sum(constinflong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    constinfshort$dose4_selec[k]<-sum(constinflong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    constinfshort$dose5_selec[k]<-sum(constinflong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    constinfshort$dose1_treat[k]<- sum(constinflong$numdose1[((k-1)*100+1):(k*100)])*100/sum(constinflong$n[((k-1)*100+1):(k*100)])
    constinfshort$dose2_treat[k]<- sum(constinflong$numdose2[((k-1)*100+1):(k*100)])*100/sum(constinflong$n[((k-1)*100+1):(k*100)])
    constinfshort$dose3_treat[k]<- sum(constinflong$numdose3[((k-1)*100+1):(k*100)])*100/sum(constinflong$n[((k-1)*100+1):(k*100)])
    constinfshort$dose4_treat[k]<- sum(constinflong$numdose4[((k-1)*100+1):(k*100)])*100/sum(constinflong$n[((k-1)*100+1):(k*100)])
    constinfshort$dose5_treat[k]<- sum(constinflong$numdose5[((k-1)*100+1):(k*100)])*100/sum(constinflong$n[((k-1)*100+1):(k*100)])
    
    constinfshort$avgDLT[k] <- mean(constinflong$numDLT[((k-1)*100+1):(k*100)])
    constinfshort$avgDR[k] <- mean(constinflong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(constinfshort,paste("constinf_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for constant
c1 <- read.csv("constinf_summaryresults_1.csv")
c2<- read.csv("constinf_summaryresults_2.csv")
c3<- read.csv("constinf_summaryresults_3.csv")
c4<- read.csv("constinf_summaryresults_4.csv")
c5<- read.csv("constinf_summaryresults_5.csv")
c6 <- read.csv("constinf_summaryresults_6.csv")
c7<- read.csv("constinf_summaryresults_7.csv")
c8<- read.csv("constinf_summaryresults_8.csv")
c9<- read.csv("constinf_summaryresults_9.csv")
c10<- read.csv("constinf_summaryresults_10.csv")


constinfshort <- data.frame(matrix(nrow=54,ncol=16))
names(constinfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                          "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
constinfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
constinfshort$targetdose <- rep(c(1,3,5),18)
constinfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

constinfshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10
constinfshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                            c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

constinfshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                           c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

constinfshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                              c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
constinfshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                              c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
constinfshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                              c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
constinfshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                              c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
constinfshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                              c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

constinfshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                               c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
constinfshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                               c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
constinfshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                               c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
constinfshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                               c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
constinfshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                               c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10

write.csv(constinfshort,file="constinf_summaryresults.csv")








###### Uninform Normal prior

####################Constant generation, Constant analysis. 

for(file in 1:10){
  load(paste("results_constgen_constanalysis_uninform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  const_uninflong <- data.frame(matrix(nrow=5400,ncol=12))
  names(const_uninflong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  const_uninflong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  const_uninflong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  const_uninflong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      const_uninflong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      const_uninflong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      const_uninflong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      const_uninflong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      const_uninflong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      const_uninflong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      const_uninflong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      const_uninflong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      const_uninflong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  const_uninflong$chosendoseorder <- as.integer(as.factor(const_uninflong$chosendose))
  
  
  
  const_uninfshort <- data.frame(matrix(nrow=54,ncol=16))
  names(const_uninfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                            "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  const_uninfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  const_uninfshort$targetdose <- rep(c(1,3,5),18)
  const_uninfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- const_uninflong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(const_uninflong$duration[((k-1)*100+1):(k*100)],lag=1)
    const_uninfshort$avgduration[k] <- mean(durationdiff)
    
    const_uninfshort$dose1_selec[k]<-sum(const_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    const_uninfshort$dose2_selec[k]<-sum(const_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    const_uninfshort$dose3_selec[k]<-sum(const_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    const_uninfshort$dose4_selec[k]<-sum(const_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    const_uninfshort$dose5_selec[k]<-sum(const_uninflong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    const_uninfshort$dose1_treat[k]<- sum(const_uninflong$numdose1[((k-1)*100+1):(k*100)])*100/sum(const_uninflong$n[((k-1)*100+1):(k*100)])
    const_uninfshort$dose2_treat[k]<- sum(const_uninflong$numdose2[((k-1)*100+1):(k*100)])*100/sum(const_uninflong$n[((k-1)*100+1):(k*100)])
    const_uninfshort$dose3_treat[k]<- sum(const_uninflong$numdose3[((k-1)*100+1):(k*100)])*100/sum(const_uninflong$n[((k-1)*100+1):(k*100)])
    const_uninfshort$dose4_treat[k]<- sum(const_uninflong$numdose4[((k-1)*100+1):(k*100)])*100/sum(const_uninflong$n[((k-1)*100+1):(k*100)])
    const_uninfshort$dose5_treat[k]<- sum(const_uninflong$numdose5[((k-1)*100+1):(k*100)])*100/sum(const_uninflong$n[((k-1)*100+1):(k*100)])
    
    const_uninfshort$avgDLT[k] <- mean(const_uninflong$numDLT[((k-1)*100+1):(k*100)])
    const_uninfshort$avgDR[k] <- mean(const_uninflong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(const_uninfshort,paste("const_uninf_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for costant
c1 <- read.csv("const_uninf_summaryresults_1.csv")
c2<- read.csv("const_uninf_summaryresults_2.csv")
c3<- read.csv("const_uninf_summaryresults_3.csv")
c4<- read.csv("const_uninf_summaryresults_4.csv")
c5<- read.csv("const_uninf_summaryresults_5.csv")
c6 <- read.csv("const_uninf_summaryresults_6.csv")
c7<- read.csv("const_uninf_summaryresults_7.csv")
c8<- read.csv("const_uninf_summaryresults_8.csv")
c9<- read.csv("const_uninf_summaryresults_9.csv")
c10<- read.csv("const_uninf_summaryresults_10.csv")


const_uninfshort <- data.frame(matrix(nrow=54,ncol=16))
names(const_uninfshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                          "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
const_uninfshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
const_uninfshort$targetdose <- rep(c(1,3,5),18)
const_uninfshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

const_uninfshort$avgduration <- (c1$avgduration+c2$avgduration+c3$avgduration+c4$avgduration+c5$avgduration+
                                c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/10
const_uninfshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c3$avgDLT+c4$avgDLT+c5$avgDLT+
                           c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/10

const_uninfshort$avgDR <- (c1$avgDR+c2$avgDR+c3$avgDR+c4$avgDR+c5$avgDR+
                          c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/10

const_uninfshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c3$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                              c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/10
const_uninfshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c3$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                              c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/10
const_uninfshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c3$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                              c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/10
const_uninfshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c3$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                              c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/10
const_uninfshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c3$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                              c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/10

const_uninfshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c3$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                               c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/10
const_uninfshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c3$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                               c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/10
const_uninfshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c3$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                               c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/10
const_uninfshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c3$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                               c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/10
const_uninfshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c3$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                               c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/10

write.csv(const_uninfshort,file="const_uninf_summaryresults.csv")














###### Uniform prior

####################Constant generation, Constant analysis. 


for(file in c(1,2,4:10)){
  load(paste("results_constgen_constanalysis_uniformprior_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  const_uniformpriorlong <- data.frame(matrix(nrow=5400,ncol=12))
  names(const_uniformpriorlong) <- c("dlttarget","targetdose","chosendose","duration","n","numDLT","numDR","numdose1","numdose2","numdose3","numdose4","numdose5")
  
  const_uniformpriorlong$dlttarget <- rep(c(rep(30,300),rep(20,300)),9)
  const_uniformpriorlong$targetdose <- rep(c(rep(1,100),rep(3,100),rep(5,100)),18)
  const_uniformpriorlong$n <- c(rep(100,1800),rep(60,1800),rep(30,1800))
  
  
  for(i in 1:54){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      const_uniformpriorlong$chosendose[(i-1)*100+j] <- spec$dose0[nrow(spec)]
      const_uniformpriorlong$duration[(i-1)*100+j] <- spec$accrualtime[nrow(spec)]
      const_uniformpriorlong$numdose1[(i-1)*100+j] <- sum(spec$dose0==0.1)
      const_uniformpriorlong$numdose2[(i-1)*100+j] <- sum(spec$dose0==0.2)
      const_uniformpriorlong$numdose3[(i-1)*100+j] <- sum(spec$dose0==0.4)
      const_uniformpriorlong$numdose4[(i-1)*100+j] <- sum(spec$dose0==0.65)
      const_uniformpriorlong$numdose5[(i-1)*100+j] <- sum(spec$dose0==1)
      const_uniformpriorlong$numDLT[(i-1)*100+j] <- spec$total_dlt[nrow(spec)]
      const_uniformpriorlong$numDR[(i-1)*100+j] <- spec$total_dr[nrow(spec)]
    }
  }
  const_uniformpriorlong$chosendoseorder <- as.integer(as.factor(const_uniformpriorlong$chosendose))
  
  
  
  const_uniformpriorshort <- data.frame(matrix(nrow=54,ncol=16))
  names(const_uniformpriorshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                               "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
  const_uniformpriorshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
  const_uniformpriorshort$targetdose <- rep(c(1,3,5),18)
  const_uniformpriorshort$n <- c(rep(100,18),rep(60,18),rep(30,18))
  
  durationdiff <- vector()
  for(k in 1:54){
    durationdiff[1]<- const_uniformpriorlong$duration[((k-1)*100+1)]
    durationdiff[2:100] <- diff(const_uniformpriorlong$duration[((k-1)*100+1):(k*100)],lag=1)
    const_uniformpriorshort$avgduration[k] <- mean(durationdiff)
    
    const_uniformpriorshort$dose1_selec[k]<-sum(const_uniformpriorlong$chosendoseorder[((k-1)*100+1):(k*100)]==1)*100/100
    const_uniformpriorshort$dose2_selec[k]<-sum(const_uniformpriorlong$chosendoseorder[((k-1)*100+1):(k*100)]==2)*100/100
    const_uniformpriorshort$dose3_selec[k]<-sum(const_uniformpriorlong$chosendoseorder[((k-1)*100+1):(k*100)]==3)*100/100
    const_uniformpriorshort$dose4_selec[k]<-sum(const_uniformpriorlong$chosendoseorder[((k-1)*100+1):(k*100)]==4)*100/100
    const_uniformpriorshort$dose5_selec[k]<-sum(const_uniformpriorlong$chosendoseorder[((k-1)*100+1):(k*100)]==5)*100/100
    
    const_uniformpriorshort$dose1_treat[k]<- sum(const_uniformpriorlong$numdose1[((k-1)*100+1):(k*100)])*100/sum(const_uniformpriorlong$n[((k-1)*100+1):(k*100)])
    const_uniformpriorshort$dose2_treat[k]<- sum(const_uniformpriorlong$numdose2[((k-1)*100+1):(k*100)])*100/sum(const_uniformpriorlong$n[((k-1)*100+1):(k*100)])
    const_uniformpriorshort$dose3_treat[k]<- sum(const_uniformpriorlong$numdose3[((k-1)*100+1):(k*100)])*100/sum(const_uniformpriorlong$n[((k-1)*100+1):(k*100)])
    const_uniformpriorshort$dose4_treat[k]<- sum(const_uniformpriorlong$numdose4[((k-1)*100+1):(k*100)])*100/sum(const_uniformpriorlong$n[((k-1)*100+1):(k*100)])
    const_uniformpriorshort$dose5_treat[k]<- sum(const_uniformpriorlong$numdose5[((k-1)*100+1):(k*100)])*100/sum(const_uniformpriorlong$n[((k-1)*100+1):(k*100)])
    
    const_uniformpriorshort$avgDLT[k] <- mean(const_uniformpriorlong$numDLT[((k-1)*100+1):(k*100)])
    const_uniformpriorshort$avgDR[k] <- mean(const_uniformpriorlong$numDR[((k-1)*100+1):(k*100)])
  }
  
  
  write.csv(const_uniformpriorshort,paste("const_uniformprior_summaryresults_",file,".csv",sep=''))
}


###### Put together final averages for costant
c1 <- read.csv("const_uniformprior_summaryresults_1.csv")
c2<- read.csv("const_uniformprior_summaryresults_2.csv")
#c3<- read.csv("const_uniformprior_summaryresults_3.csv")
c4<- read.csv("const_uniformprior_summaryresults_4.csv")
c5<- read.csv("const_uniformprior_summaryresults_5.csv")
c6 <- read.csv("const_uniformprior_summaryresults_6.csv")
c7<- read.csv("const_uniformprior_summaryresults_7.csv")
c8<- read.csv("const_uniformprior_summaryresults_8.csv")
c9<- read.csv("const_uniformprior_summaryresults_9.csv")
c10<- read.csv("const_uniformprior_summaryresults_10.csv")


const_uniformpriorshort <- data.frame(matrix(nrow=54,ncol=16))
names(const_uniformpriorshort) <- c("dlttarget","targetdose","avgduration","n","avgDLT","avgDR","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec",
                             "dose1_treat","dose2_treat","dose3_treat","dose4_treat","dose5_treat")
const_uniformpriorshort$dlttarget <- rep(c(rep(30,3),rep(20,3)),9)
const_uniformpriorshort$targetdose <- rep(c(1,3,5),18)
const_uniformpriorshort$n <- c(rep(100,18),rep(60,18),rep(30,18))

const_uniformpriorshort$avgduration <- (c1$avgduration+c2$avgduration+c4$avgduration+c5$avgduration+
                                   c6$avgduration+c7$avgduration+c8$avgduration+c9$avgduration+c10$avgduration)/9
const_uniformpriorshort$avgDLT <- (c1$avgDLT+c2$avgDLT+c4$avgDLT+c5$avgDLT+
                              c6$avgDLT+c7$avgDLT+c8$avgDLT+c9$avgDLT+c10$avgDLT)/9

const_uniformpriorshort$avgDR <- (c1$avgDR+c2$avgDR+c4$avgDR+c5$avgDR+
                             c6$avgDR+c7$avgDR+c8$avgDR+c9$avgDR+c10$avgDR)/9

const_uniformpriorshort$dose1_selec<-(c1$dose1_selec+c2$dose1_selec+c4$dose1_selec+c5$dose1_selec+
                                 c6$dose1_selec+c7$dose1_selec+c8$dose1_selec+c9$dose1_selec+c10$dose1_selec)/9
const_uniformpriorshort$dose2_selec<-(c1$dose2_selec+c2$dose2_selec+c4$dose2_selec+c5$dose2_selec+
                                 c6$dose2_selec+c7$dose2_selec+c8$dose2_selec+c9$dose2_selec+c10$dose2_selec)/9
const_uniformpriorshort$dose3_selec<-(c1$dose3_selec+c2$dose3_selec+c4$dose3_selec+c5$dose3_selec+
                                 c6$dose3_selec+c7$dose3_selec+c8$dose3_selec+c9$dose3_selec+c10$dose3_selec)/9
const_uniformpriorshort$dose4_selec<-(c1$dose4_selec+c2$dose4_selec+c4$dose4_selec+c5$dose4_selec+
                                 c6$dose4_selec+c7$dose4_selec+c8$dose4_selec+c9$dose4_selec+c10$dose4_selec)/9
const_uniformpriorshort$dose5_selec<-(c1$dose5_selec+c2$dose5_selec+c4$dose5_selec+c5$dose5_selec+
                                 c6$dose5_selec+c7$dose5_selec+c8$dose5_selec+c9$dose5_selec+c10$dose5_selec)/9

const_uniformpriorshort$dose1_treat<- (c1$dose1_treat+c2$dose1_treat+c4$dose1_treat+c5$dose1_treat+
                                  c6$dose1_treat+c7$dose1_treat+c8$dose1_treat+c9$dose1_treat+c10$dose1_treat)/9
const_uniformpriorshort$dose2_treat<- (c1$dose2_treat+c2$dose2_treat+c4$dose2_treat+c5$dose2_treat+
                                  c6$dose2_treat+c7$dose2_treat+c8$dose2_treat+c9$dose2_treat+c10$dose2_treat)/9
const_uniformpriorshort$dose3_treat<- (c1$dose3_treat+c2$dose3_treat+c4$dose3_treat+c5$dose3_treat+
                                  c6$dose3_treat+c7$dose3_treat+c8$dose3_treat+c9$dose3_treat+c10$dose3_treat)/9
const_uniformpriorshort$dose4_treat<- (c1$dose4_treat+c2$dose4_treat+c4$dose4_treat+c5$dose4_treat+
                                  c6$dose4_treat+c7$dose4_treat+c8$dose4_treat+c9$dose4_treat+c10$dose4_treat)/9
const_uniformpriorshort$dose5_treat<- (c1$dose5_treat+c2$dose5_treat+c4$dose5_treat+c5$dose5_treat+
                                  c6$dose5_treat+c7$dose5_treat+c8$dose5_treat+c9$dose5_treat+c10$dose5_treat)/9

write.csv(const_uniformpriorshort,file="const_uniformprior_summaryresults.csv")
