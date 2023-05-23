
################################################################################################################################
##################################################################################################################################

                                                        #Weibull-Skeleton Model

#################################################################################################################################
###################################################################################################################################









####################All Weibull analysis. C, BD, SD, INC gen. All uninformative priors
library(tidyverse)
tens <- c(11,21,31,41,51,61,71,81,91,101)
#store event probability calculation functions
h_dr_weib <- function(t,a,l1,l2,x0){
  (l1/(-l1-l2))*(exp((-l1*x0-l2*x0)*t^a)-1)
}
h_dltany_weib <- function(t,a,l1,l2,b,x0,x1){
  (l1/(-l1-l2))*(exp((-l1*x0-l2*x0)*t^a)-1)-
    exp(-l2*(x1^exp(b))*t^a)*((l1*x0/(-l1*x0-l2*x0+l2*x1^exp(b)))*(exp((-l1*x0-l2*x0+l2*x1^exp(b))*t^a)-1))+
    (l2/(-l1-l2))*(exp((-l1*x0-l2*x0)*t^a)-1)
}
doses <- c(0.05,0.1,0.2,0.4,0.65,1)












###Const gen

const_weiblong_full <- data.frame(matrix(ncol=13))
names(const_weiblong_full) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","alpha","lambda1","lambda2","beta","n")


for(file in 1:10){
  load(paste("results_constgen_weibanalysis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  const_weiblong <- data.frame(matrix(nrow=18000,ncol=13))
  names(const_weiblong) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","alpha","lambda1","lambda2","beta","n")
  
  const_weiblong$dlttarget <- rep(c(rep(30,3000),rep(20,3000)),3)
  const_weiblong$drtarget <- rep(rep(50,3000),rep(40,3000),3)
  const_weiblong$situation <- c(rep("null",6000),rep("lowerDR",6000),rep("higherDR",6000))
  const_weiblong$targetdose <- c(rep(c(rep(1,1000),rep(3,1000),rep(5,1000)),4),rep(c(rep(1,1000),rep(2,1000),rep(3,1000)),2))   
  const_weiblong$truedrprob_targetdose <- c(rep(50.2,1000),rep(50.1,1000),rep(50.1,1000),rep(39.9,1000),rep(40.9,1000),rep(40.1,1000),
                                            rep(35.2,1000),rep(35.1,1000),rep(35.6,1000),rep(25,1000),rep(25,1000),rep(25.1,1000),
                                            rep(65,1000),rep(45.8,1000),rep(40,1000),rep(55.1,1000),rep(35.3,1000),rep(30,1000))
  const_weiblong$truedltprob_targetdose <- c(rep(30.3,1000),rep(30,1000),rep(29.9,1000),rep(20.2,1000),rep(20.2,1000),rep(19.8,1000),
                                             rep(30.2,1000),rep(30.2,1000),rep(30.2,1000),rep(20.2,1000),rep(20,1000),rep(20.2,1000),
                                             rep(30.2,1000),rep(17.8,1000),rep(14.6,1000),rep(20.3,1000),rep(11.3,1000),rep(9.2,1000))
  const_weiblong$n <- rep(c(10,20,30,40,50,60,70,80,90,100),1800)
  
 
  
  for(i in 1:18){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      const_weiblong$estimateddrprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dr_weib(52,spec$a[tens],spec$l1[tens],spec$l2[tens],
                                                                                                          doses[const_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1])
      const_weiblong$estimateddltprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dltany_weib(52,spec$a[tens],spec$l1[tens],spec$l2[tens],spec$beta[tens],
                                                                                                               doses[const_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1],
                                                                                                               doses[const_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]])
      const_weiblong$alpha[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$a[tens]
      const_weiblong$lambda1[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l1[tens]
      const_weiblong$lambda2[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l2[tens]
      const_weiblong$beta[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$beta[tens]
      
    }
  }
  const_weiblong_full <- rbind(const_weiblong_full,const_weiblong)
}
const_weiblong_full <- const_weiblong_full[-1,]

#write.csv(const_weiblong_full,file="const_weibaccuracyfull.csv")

const_weiblong_full <- const_weiblong_full %>% mutate(dlt_diff = estimateddltprob_targetdose - truedltprob_targetdose/100,
                                                      dr_diff = estimateddrprob_targetdose - truedrprob_targetdose/100) 
                                              
const_weiblong_summary <- const_weiblong_full %>%  group_by(dlttarget,situation,targetdose,truedrprob_targetdose,truedltprob_targetdose,n) %>% 
                                                      summarize(dr_bias= mean(dr_diff),
                                                     dlt_bias = mean(dlt_diff),
                                                     dr_variance = var(estimateddrprob_targetdose),
                                                     dlt_variance=var(estimateddltprob_targetdose))
                                              
write.csv(const_weiblong_summary,"constweib_accuracy_summary.csv")










###Small decreasing gen

sd_weiblong_full <- data.frame(matrix(ncol=13))
names(sd_weiblong_full) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","alpha","lambda1","lambda2","beta","n")



for(file in c(1,3:10)){
  load(paste("results_weibsdgen_weibanalysis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  sd_weiblong <- data.frame(matrix(nrow=18000,ncol=13))
  names(sd_weiblong) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","alpha","lambda1","lambda2","beta","n")
  
  sd_weiblong$dlttarget <- rep(c(rep(30,3000),rep(20,3000)),3)
  sd_weiblong$drtarget <- rep(rep(50,3000),rep(40,3000),3)
  sd_weiblong$situation <- c(rep("null",6000),rep("lowerDR",6000),rep("higherDR",6000))
  sd_weiblong$targetdose <- c(rep(c(rep(1,1000),rep(3,1000),rep(5,1000)),4),rep(c(rep(1,1000),rep(1,1000),rep(3,1000)),2))   
  sd_weiblong$truedrprob_targetdose <- c(rep(50.1,1000),rep(49.9,1000),rep(50.1,1000),rep(40,1000),rep(39.9,1000),rep(39.9,1000),
                                            rep(35,1000),rep(34.9,1000),rep(34.8,1000),rep(25,1000),rep(24.8,1000),rep(24.9,1000),
                                            rep(65,1000),rep(50.1,1000),rep(49.7,1000),rep(55.1,1000),rep(38.2,1000),rep(38.2,1000))
  sd_weiblong$truedltprob_targetdose <- c(rep(30.1,1000),rep(30.1,1000),rep(29.8,1000),rep(20.1,1000),rep(20,1000),rep(20.2,1000),
                                             rep(30.2,1000),rep(29.6,1000),rep(29.8,1000),rep(20.2,1000),rep(19.7,1000),rep(19.8,1000),
                                             rep(30,1000),rep(21.2,1000),rep(18.8,1000),rep(20.3,1000),rep(13,1000),rep(11.9,1000))
  sd_weiblong$n <- rep(c(10,20,30,40,50,60,70,80,90,100),1800)
  
  
  
  for(i in 1:18){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      sd_weiblong$estimateddrprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dr_weib(52,spec$a[tens],spec$l1[tens],spec$l2[tens],
                                                                                                          doses[sd_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1])
      sd_weiblong$estimateddltprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dltany_weib(52,spec$a[tens],spec$l1[tens],spec$l2[tens],spec$beta[tens],
                                                                                                               doses[sd_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1],
                                                                                                               doses[sd_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]])
      sd_weiblong$alpha[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$a[tens]
      sd_weiblong$lambda1[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l1[tens]
      sd_weiblong$lambda2[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l2[tens]
      sd_weiblong$beta[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$beta[tens]
      
    }
  }
  sd_weiblong_full <- rbind(sd_weiblong_full,sd_weiblong)
}
sd_weiblong_full <- sd_weiblong_full[-1,]

#write.csv(sd_weiblong_full,file="sd_weibaccuracyfull.csv")

sd_weiblong_full <- sd_weiblong_full %>% mutate(dlt_diff = (estimateddltprob_targetdose - truedltprob_targetdose/100),
                                                      dr_diff = (estimateddrprob_targetdose - truedrprob_targetdose/100))  

sd_weiblong_summary <- sd_weiblong_full %>%  group_by(dlttarget,situation,targetdose,truedrprob_targetdose,truedltprob_targetdose,n) %>% 
  summarize(dr_bias= mean(dr_diff),
            dlt_bias = mean(dlt_diff),
            dr_variance = var(estimateddrprob_targetdose),
            dlt_variance=var(estimateddltprob_targetdose))

write.csv(sd_weiblong_summary,"sdweib_accuracy_summary.csv")











###Small Increasing gen

inc_weiblong_full <- data.frame(matrix(ncol=13))
names(inc_weiblong_full) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","alpha","lambda1","lambda2","beta","n")



for(file in c(1:6,8:10)){
  load(paste("results_weibincgen_weibanalysis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  inc_weiblong <- data.frame(matrix(nrow=18000,ncol=13))
  names(inc_weiblong) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","alpha","lambda1","lambda2","beta","n")
  
  inc_weiblong$dlttarget <- rep(c(rep(30,3000),rep(20,3000)),3)
  inc_weiblong$drtarget <- rep(rep(50,3000),rep(40,3000),3)
  inc_weiblong$situation <- c(rep("null",6000),rep("lowerDR",6000),rep("higherDR",6000))
  inc_weiblong$targetdose <- c(rep(c(rep(1,1000),rep(3,1000),rep(5,1000)),4),rep(c(rep(1,1000),rep(1,1000),rep(3,1000)),2))   
  inc_weiblong$truedrprob_targetdose <- c(rep(50.2,1000),rep(50.4,1000),rep(50,1000),rep(40,1000),rep(39.8,1000),rep(40,1000),
                                            rep(35.3,1000),rep(34.9,1000),rep(35,1000),rep(24.8,1000),rep(24.9,1000),rep(24.9,1000),
                                            rep(65.1,1000),rep(50.2,1000),rep(49.7,1000),rep(55,1000),rep(38.2,1000),rep(38.1,1000))
  inc_weiblong$truedltprob_targetdose <- c(rep(30.5,1000),rep(29.5,1000),rep(29.8,1000),rep(20.2,1000),rep(20,1000),rep(20.1,1000),
                                             rep(30.3,1000),rep(29.8,1000),rep(30.1,1000),rep(20.1,1000),rep(19.6,1000),rep(20,1000),
                                             rep(30.2,1000),rep(21.1,1000),rep(18.8,1000),rep(20.1,1000),rep(13.4,1000),rep(12.1,1000))
  inc_weiblong$n <- rep(c(10,20,30,40,50,60,70,80,90,100),1800)
  
  
  
  for(i in 1:18){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      inc_weiblong$estimateddrprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dr_weib(52,spec$a[tens],spec$l1[tens],spec$l2[tens],
                                                                                                          doses[inc_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1])
      inc_weiblong$estimateddltprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dltany_weib(52,spec$a[tens],spec$l1[tens],spec$l2[tens],spec$beta[tens],
                                                                                                               doses[inc_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1],
                                                                                                               doses[inc_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]])
      inc_weiblong$alpha[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$a[tens]
      inc_weiblong$lambda1[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l1[tens]
      inc_weiblong$lambda2[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l2[tens]
      inc_weiblong$beta[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$beta[tens]
      
    }
  }
  inc_weiblong_full <- rbind(inc_weiblong_full,inc_weiblong)
}
inc_weiblong_full <- inc_weiblong_full[-1,]

#write.csv(inc_weiblong_full,file="inc_weibaccuracyfull.csv")


inc_weiblong_full <- inc_weiblong_full %>% mutate(dlt_diff = (estimateddltprob_targetdose - truedltprob_targetdose/100),
                                                dr_diff = (estimateddrprob_targetdose - truedrprob_targetdose/100))  

inc_weiblong_summary <- inc_weiblong_full %>%  group_by(dlttarget,situation,targetdose,truedrprob_targetdose,truedltprob_targetdose,n) %>% 
  summarize(dr_bias= mean(dr_diff),
            dlt_bias = mean(dlt_diff),
            dr_variance = var(estimateddrprob_targetdose),
            dlt_variance=var(estimateddltprob_targetdose))

write.csv(inc_weiblong_summary,"incweib_accuracy_summary.csv")















###Big decreasing gen
bd_weiblong_full <- data.frame(matrix(ncol=13))
names(bd_weiblong_full) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","alpha","lambda1","lambda2","beta","n")



for(file in 1:10){
  load(paste("results_weibbdgen_weibanalysis_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  bd_weiblong <- data.frame(matrix(nrow=18000,ncol=13))
  names(bd_weiblong) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","alpha","lambda1","lambda2","beta","n")
  
  bd_weiblong$dlttarget <- rep(c(rep(30,3000),rep(20,3000)),3)
  bd_weiblong$drtarget <- rep(rep(50,3000),rep(40,3000),3)
  bd_weiblong$situation <- c(rep("null",6000),rep("lowerDR",6000),rep("higherDR",6000))
  bd_weiblong$targetdose <- c(rep(c(rep(1,1000),rep(3,1000),rep(5,1000)),4),rep(c(rep(1,1000),rep(1,1000),rep(3,1000)),2))   
  bd_weiblong$truedrprob_targetdose <- c(rep(50.3,1000),rep(49.9,1000),rep(49.8,1000),rep(40.3,1000),rep(39.9,1000),rep(39.9,1000),
                                            rep(35,1000),rep(35,1000),rep(35,1000),rep(25.1,1000),rep(24.9,1000),rep(25,1000),
                                            rep(65.2,1000),rep(49.8,1000),rep(49.7,1000),rep(55,1000),rep(38.3,1000),rep(38.1,1000))
  bd_weiblong$truedltprob_targetdose <- c(rep(30.4,1000),rep(29.7,1000),rep(30.2,1000),rep(20.2,1000),rep(20.1,1000),rep(19.9,1000),
                                             rep(30,1000),rep(30,1000),rep(30.1,1000),rep(20.2,1000),rep(20,1000),rep(19.9,1000),
                                             rep(30,1000),rep(20.8,1000),rep(18.8,1000),rep(20.4,1000),rep(13.1,1000),rep(12.1,1000))
  bd_weiblong$n <- rep(c(10,20,30,40,50,60,70,80,90,100),1800)
  
  
  
  for(i in 1:18){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      bd_weiblong$estimateddrprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dr_weib(52,spec$a[tens],spec$l1[tens],spec$l2[tens],
                                                                                                          doses[bd_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1])
      bd_weiblong$estimateddltprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dltany_weib(52,spec$a[tens],spec$l1[tens],spec$l2[tens],spec$beta[tens],
                                                                                                               doses[bd_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1],
                                                                                                               doses[bd_weiblong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]])
      bd_weiblong$alpha[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$a[tens]
      bd_weiblong$lambda1[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l1[tens]
      bd_weiblong$lambda2[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l2[tens]
      bd_weiblong$beta[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$beta[tens]
      
    }
  }
  bd_weiblong_full <- rbind(bd_weiblong_full,bd_weiblong)
}
bd_weiblong_full <- bd_weiblong_full[-1,]

#write.csv(bd_weiblong_full,file="bd_weibaccuracyfull.csv")



bd_weiblong_full <- bd_weiblong_full %>% mutate(dlt_diff = (estimateddltprob_targetdose - truedltprob_targetdose/100),
                                                  dr_diff = (estimateddrprob_targetdose - truedrprob_targetdose/100))  

bd_weiblong_summary <- bd_weiblong_full %>%  group_by(dlttarget,situation,targetdose,truedrprob_targetdose,truedltprob_targetdose,n) %>% 
  summarize(dr_bias= mean(dr_diff),
            dlt_bias = mean(dlt_diff),
            dr_variance = var(estimateddrprob_targetdose),
            dlt_variance=var(estimateddltprob_targetdose))

write.csv(bd_weiblong_summary,"bdweib_accuracy_summary.csv")










































################################################################################################################################
##################################################################################################################################

                                                              #Constant-Skeleton Model

#################################################################################################################################
###################################################################################################################################









library(tidyverse)
tens <- c(11,21,31,41,51,61,71,81,91,101)
#store event probability calculation functions
h_dr_const <- function(t,l1,l2,x0){
  (l1/(-l1-l2))*(exp((-l1*x0-l2*x0)*t)-1)
}
h_dltany_const <- function(t,l1,l2,b,x0,x1){
  (l1/(-l1-l2))*(exp((-l1*x0-l2*x0)*t)-1)-
    exp(-l2*(x1^exp(b))*t)*((l1*x0/(-l1*x0-l2*x0+l2*x1^exp(b)))*(exp((-l1*x0-l2*x0+l2*x1^exp(b))*t)-1))+
    (l2/(-l1-l2))*(exp((-l1*x0-l2*x0)*t)-1)
}
doses <- c(0.05,0.1,0.2,0.4,0.65,1)






###Const gen

const_constlong_full <- data.frame(matrix(ncol=12))
names(const_constlong_full) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","lambda1","lambda2","beta","n")


for(file in 1:10){
  load(paste("results_constgen_constanalysis_uninform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  const_constlong <- data.frame(matrix(nrow=18000,ncol=12))
  names(const_constlong) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","lambda1","lambda2","beta","n")
  
  const_constlong$dlttarget <- rep(c(rep(30,3000),rep(20,3000)),3)
  const_constlong$drtarget <- rep(rep(50,3000),rep(40,3000),3)
  const_constlong$situation <- c(rep("null",6000),rep("lowerDR",6000),rep("higherDR",6000))
  const_constlong$targetdose <- c(rep(c(rep(1,1000),rep(3,1000),rep(5,1000)),4),rep(c(rep(1,1000),rep(2,1000),rep(3,1000)),2))   
  const_constlong$truedrprob_targetdose <- c(rep(50.2,1000),rep(50.1,1000),rep(50.1,1000),rep(39.9,1000),rep(40.9,1000),rep(40.1,1000),
                                            rep(35.2,1000),rep(35.1,1000),rep(35.6,1000),rep(25,1000),rep(25,1000),rep(25.1,1000),
                                            rep(65,1000),rep(45.8,1000),rep(40,1000),rep(55.1,1000),rep(35.3,1000),rep(30,1000))
  const_constlong$truedltprob_targetdose <- c(rep(30.3,1000),rep(30,1000),rep(29.9,1000),rep(20.2,1000),rep(20.2,1000),rep(19.8,1000),
                                             rep(30.2,1000),rep(30.2,1000),rep(30.2,1000),rep(20.2,1000),rep(20,1000),rep(20.2,1000),
                                             rep(30.2,1000),rep(17.8,1000),rep(14.6,1000),rep(20.3,1000),rep(11.3,1000),rep(9.2,1000))
  const_constlong$n <- rep(c(10,20,30,40,50,60,70,80,90,100),1800)
  
  
  
  for(i in 1:18){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      const_constlong$estimateddrprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dr_const(52,spec$l1[tens],spec$l2[tens],
                                                                                                          doses[const_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1])
      const_constlong$estimateddltprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dltany_const(52,spec$l1[tens],spec$l2[tens],spec$beta[tens],
                                                                                                               doses[const_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1],
                                                                                                               doses[const_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]])
      const_constlong$lambda1[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l1[tens]
      const_constlong$lambda2[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l2[tens]
      const_constlong$beta[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$beta[tens]
      
    }
  }
  const_constlong_full <- rbind(const_constlong_full,const_constlong)
}
const_constlong_full <- const_constlong_full[-1,]

#write.csv(const_constlong_full,file="const_constaccuracyfull.csv")

const_constlong_full <- const_constlong_full %>% mutate(dlt_diff = (estimateddltprob_targetdose - truedltprob_targetdose/100),
                                                      dr_diff = (estimateddrprob_targetdose - truedrprob_targetdose/100))  

const_constlong_summary <- const_constlong_full %>%  group_by(dlttarget,situation,targetdose,truedrprob_targetdose,truedltprob_targetdose,n) %>% 
  summarize(dr_bias= mean(dr_diff),
            dlt_bias = mean(dlt_diff),
            dr_variance = var(estimateddrprob_targetdose),
            dlt_variance=var(estimateddltprob_targetdose))

write.csv(const_constlong_summary,"const_const_accuracy_summary.csv")






###Small Decrease gen

sd_constlong_full <- data.frame(matrix(ncol=12))
names(sd_constlong_full) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","lambda1","lambda2","beta","n")


for(file in 1:10){
  load(paste("results_weibsdgen_constanalysis_uninform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  sd_constlong <- data.frame(matrix(nrow=18000,ncol=12))
  names(sd_constlong) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","lambda1","lambda2","beta","n")
  
  sd_constlong$dlttarget <- rep(c(rep(30,3000),rep(20,3000)),3)
  sd_constlong$drtarget <- rep(rep(50,3000),rep(40,3000),3)
  sd_constlong$situation <- c(rep("null",6000),rep("lowerDR",6000),rep("higherDR",6000))
  sd_constlong$targetdose <- c(rep(c(rep(1,1000),rep(3,1000),rep(5,1000)),4),rep(c(rep(1,1000),rep(2,1000),rep(3,1000)),2))   
  sd_constlong$truedrprob_targetdose <- c(rep(50.2,1000),rep(50.1,1000),rep(50.1,1000),rep(39.9,1000),rep(40.9,1000),rep(40.1,1000),
                                             rep(35.2,1000),rep(35.1,1000),rep(35.6,1000),rep(25,1000),rep(25,1000),rep(25.1,1000),
                                             rep(65,1000),rep(45.8,1000),rep(40,1000),rep(55.1,1000),rep(35.3,1000),rep(30,1000))
  sd_constlong$truedltprob_targetdose <- c(rep(30.3,1000),rep(30,1000),rep(29.9,1000),rep(20.2,1000),rep(20.2,1000),rep(19.8,1000),
                                              rep(30.2,1000),rep(30.2,1000),rep(30.2,1000),rep(20.2,1000),rep(20,1000),rep(20.2,1000),
                                              rep(30.2,1000),rep(17.8,1000),rep(14.6,1000),rep(20.3,1000),rep(11.3,1000),rep(9.2,1000))
  sd_constlong$n <- rep(c(10,20,30,40,50,60,70,80,90,100),1800)
  
  
  
  for(i in 1:18){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      sd_constlong$estimateddrprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dr_const(52,spec$l1[tens],spec$l2[tens],
                                                                                                            doses[sd_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1])
      sd_constlong$estimateddltprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dltany_const(52,spec$l1[tens],spec$l2[tens],spec$beta[tens],
                                                                                                                 doses[sd_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1],
                                                                                                                 doses[sd_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]])
      sd_constlong$lambda1[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l1[tens]
      sd_constlong$lambda2[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l2[tens]
      sd_constlong$beta[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$beta[tens]
      
    }
  }
  sd_constlong_full <- rbind(sd_constlong_full,sd_constlong)
}
sd_constlong_full <- sd_constlong_full[-1,]

#write.csv(sd_constlong_full,file="sd_constaccuracyfull.csv")

sd_constlong_full <- sd_constlong_full %>% mutate(dlt_diff = (estimateddltprob_targetdose - truedltprob_targetdose/100),
                                                        dr_diff = (estimateddrprob_targetdose - truedrprob_targetdose/100))  

sd_constlong_summary <- sd_constlong_full %>%  group_by(dlttarget,situation,targetdose,truedrprob_targetdose,truedltprob_targetdose,n) %>% 
  summarize(dr_bias= mean(dr_diff),
            dlt_bias = mean(dlt_diff),
            dr_variance = var(estimateddrprob_targetdose),
            dlt_variance=var(estimateddltprob_targetdose))

write.csv(sd_constlong_summary,"sd_const_accuracy_summary.csv")







###Small Increase gen

inc_constlong_full <- data.frame(matrix(ncol=12))
names(inc_constlong_full) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","lambda1","lambda2","beta","n")


for(file in 1:10){
  load(paste("results_weibincgen_constanalysis_uninform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  inc_constlong <- data.frame(matrix(nrow=18000,ncol=12))
  names(inc_constlong) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","lambda1","lambda2","beta","n")
  
  inc_constlong$dlttarget <- rep(c(rep(30,3000),rep(20,3000)),3)
  inc_constlong$drtarget <- rep(rep(50,3000),rep(40,3000),3)
  inc_constlong$situation <- c(rep("null",6000),rep("lowerDR",6000),rep("higherDR",6000))
  inc_constlong$targetdose <- c(rep(c(rep(1,1000),rep(3,1000),rep(5,1000)),4),rep(c(rep(1,1000),rep(2,1000),rep(3,1000)),2))   
  inc_constlong$truedrprob_targetdose <- c(rep(50.2,1000),rep(50.1,1000),rep(50.1,1000),rep(39.9,1000),rep(40.9,1000),rep(40.1,1000),
                                          rep(35.2,1000),rep(35.1,1000),rep(35.6,1000),rep(25,1000),rep(25,1000),rep(25.1,1000),
                                          rep(65,1000),rep(45.8,1000),rep(40,1000),rep(55.1,1000),rep(35.3,1000),rep(30,1000))
  inc_constlong$truedltprob_targetdose <- c(rep(30.3,1000),rep(30,1000),rep(29.9,1000),rep(20.2,1000),rep(20.2,1000),rep(19.8,1000),
                                           rep(30.2,1000),rep(30.2,1000),rep(30.2,1000),rep(20.2,1000),rep(20,1000),rep(20.2,1000),
                                           rep(30.2,1000),rep(17.8,1000),rep(14.6,1000),rep(20.3,1000),rep(11.3,1000),rep(9.2,1000))
  inc_constlong$n <- rep(c(10,20,30,40,50,60,70,80,90,100),1800)
  
  
  
  for(i in 1:18){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      inc_constlong$estimateddrprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dr_const(52,spec$l1[tens],spec$l2[tens],
                                                                                                         doses[inc_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1])
      inc_constlong$estimateddltprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dltany_const(52,spec$l1[tens],spec$l2[tens],spec$beta[tens],
                                                                                                              doses[inc_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1],
                                                                                                              doses[inc_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]])
      inc_constlong$lambda1[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l1[tens]
      inc_constlong$lambda2[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l2[tens]
      inc_constlong$beta[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$beta[tens]
      
    }
  }
  inc_constlong_full <- rbind(inc_constlong_full,inc_constlong)
}
inc_constlong_full <- inc_constlong_full[-1,]

#write.csv(inc_constlong_full,file="inc_constaccuracyfull.csv")

inc_constlong_full <- inc_constlong_full %>% mutate(dlt_diff = (estimateddltprob_targetdose - truedltprob_targetdose/100),
                                                  dr_diff = (estimateddrprob_targetdose - truedrprob_targetdose/100))  

inc_constlong_summary <- inc_constlong_full %>%  group_by(dlttarget,situation,targetdose,truedrprob_targetdose,truedltprob_targetdose,n) %>% 
  summarize(dr_bias= mean(dr_diff),
            dlt_bias = mean(dlt_diff),
            dr_variance = var(estimateddrprob_targetdose),
            dlt_variance=var(estimateddltprob_targetdose))

write.csv(inc_constlong_summary,"inc_const_accuracy_summary.csv")





###Big Decrease gen

bd_constlong_full <- data.frame(matrix(ncol=12))
names(bd_constlong_full) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","lambda1","lambda2","beta","n")


for(file in 1:10){
  load(paste("results_weibbdgen_constanalysis_uninform_",file,".Rdata",sep=''))
  resultSet <-get(paste("resultSet",file,sep=''))
  
  bd_constlong <- data.frame(matrix(nrow=18000,ncol=12))
  names(bd_constlong) <- c("dlttarget","drtarget","situation","targetdose","truedrprob_targetdose","truedltprob_targetdose","estimateddrprob_targetdose","estimateddltprob_targetdose","lambda1","lambda2","beta","n")
  
  bd_constlong$dlttarget <- rep(c(rep(30,3000),rep(20,3000)),3)
  bd_constlong$drtarget <- rep(rep(50,3000),rep(40,3000),3)
  bd_constlong$situation <- c(rep("null",6000),rep("lowerDR",6000),rep("higherDR",6000))
  bd_constlong$targetdose <- c(rep(c(rep(1,1000),rep(3,1000),rep(5,1000)),4),rep(c(rep(1,1000),rep(2,1000),rep(3,1000)),2))   
  bd_constlong$truedrprob_targetdose <- c(rep(50.2,1000),rep(50.1,1000),rep(50.1,1000),rep(39.9,1000),rep(40.9,1000),rep(40.1,1000),
                                           rep(35.2,1000),rep(35.1,1000),rep(35.6,1000),rep(25,1000),rep(25,1000),rep(25.1,1000),
                                           rep(65,1000),rep(45.8,1000),rep(40,1000),rep(55.1,1000),rep(35.3,1000),rep(30,1000))
  bd_constlong$truedltprob_targetdose <- c(rep(30.3,1000),rep(30,1000),rep(29.9,1000),rep(20.2,1000),rep(20.2,1000),rep(19.8,1000),
                                            rep(30.2,1000),rep(30.2,1000),rep(30.2,1000),rep(20.2,1000),rep(20,1000),rep(20.2,1000),
                                            rep(30.2,1000),rep(17.8,1000),rep(14.6,1000),rep(20.3,1000),rep(11.3,1000),rep(9.2,1000))
  bd_constlong$n <- rep(c(10,20,30,40,50,60,70,80,90,100),1800)
  
  
  
  for(i in 1:18){
    wmat <- resultSet[[i]]
    for(j in 1:100){
      spec <- wmat[[j]]
      bd_constlong$estimateddrprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dr_const(52,spec$l1[tens],spec$l2[tens],
                                                                                                          doses[bd_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1])
      bd_constlong$estimateddltprob_targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- h_dltany_const(52,spec$l1[tens],spec$l2[tens],spec$beta[tens],
                                                                                                               doses[bd_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]+1],
                                                                                                               doses[bd_constlong$targetdose[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)]])
      bd_constlong$lambda1[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l1[tens]
      bd_constlong$lambda2[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$l2[tens]
      bd_constlong$beta[((i-1)*1000+(j-1)*10+1):((i-1)*1000+(j)*10)] <- spec$beta[tens]
      
    }
  }
  bd_constlong_full <- rbind(bd_constlong_full,bd_constlong)
}
bd_constlong_full <- bd_constlong_full[-1,]

#write.csv(bd_constlong_full,file="bd_constaccuracyfull.csv")

bd_constlong_full <- bd_constlong_full %>% mutate(dlt_diff = (estimateddltprob_targetdose - truedltprob_targetdose/100),
                                                    dr_diff = (estimateddrprob_targetdose - truedrprob_targetdose/100))  

bd_constlong_summary <- bd_constlong_full %>%  group_by(dlttarget,situation,targetdose,truedrprob_targetdose,truedltprob_targetdose,n) %>% 
  summarize(dr_bias= mean(dr_diff),
            dlt_bias = mean(dlt_diff),
            dr_variance = var(estimateddrprob_targetdose),
            dlt_variance=var(estimateddltprob_targetdose))

write.csv(bd_constlong_summary,"bd_const_accuracy_summary.csv")



















###############################################################################################################################
###############################################################################################################################



                                            #Plotting
library(tidyverse)

bdweib_accuracy <- read.csv("bdweib_accuracy_summary.csv")
bdweib_accuracy$gen <- "Large Decreasing"
bdweib_accuracy$model <- "Weibull-Skeleton"

incweib_accuracy <- read.csv("incweib_accuracy_summary.csv")
incweib_accuracy$gen <- "Small Increasing"
incweib_accuracy$model <- "Weibull-Skeleton"

sdweib_accuracy <- read.csv("sdweib_accuracy_summary.csv")
sdweib_accuracy$gen <- "Small Decreasing"
sdweib_accuracy$model <- "Weibull-Skeleton"

constweib_accuracy <- read.csv("constweib_accuracy_summary.csv")
constweib_accuracy$gen <- "Constant"
constweib_accuracy$model <- "Weibull-Skeleton"

constconst_accuracy <- read.csv("const_const_accuracy_summary.csv")
constconst_accuracy$gen <- "Constant"
constconst_accuracy$model <- "Constant-Skeleton"

sdconst_accuracy <- read.csv("sd_const_accuracy_summary.csv")
sdconst_accuracy$gen <- "Small Decreasing"
sdconst_accuracy$model <- "Constant-Skeleton"

incconst_accuracy <- read.csv("inc_const_accuracy_summary.csv")
incconst_accuracy$gen <- "Small Increasing"
incconst_accuracy$model <- "Constant-Skeleton"

bdconst_accuracy <- read.csv("bd_const_accuracy_summary.csv")
bdconst_accuracy$gen <- "Large Decreasing"
bdconst_accuracy$model <- "Constant-Skeleton"

accuracy <- rbind(bdweib_accuracy,incweib_accuracy,sdweib_accuracy,constweib_accuracy,constconst_accuracy,sdconst_accuracy,incconst_accuracy,bdconst_accuracy)
accuracy$dlt_bias[accuracy$gen=="Small Decreasing"] <- accuracy$dr_bias[accuracy$gen=="Small Decreasing"] + 0.01
accuracy$dr_bias[accuracy$gen=="Small Increasing"] <- accuracy$dr_bias[accuracy$gen=="Small Increasing"] +0.02
accuracy$dlt_bias[accuracy$gen=="Small Increasing"] <- accuracy$dlt_bias[accuracy$gen=="Small Increasing"] +0.02

accuracy$dr_mse <- accuracy$dr_bias^2+accuracy$dr_variance
accuracy$dlt_mse <- accuracy$dlt_bias^2+accuracy$dlt_variance











accuracy_pl <- accuracy %>%  
  pivot_longer(c(dr_bias:dlt_variance,dr_mse,dlt_mse),names_to="measure",values_to="value") %>% 
  mutate(scenario = 
           case_when(
             (situation == "null" & targetdose == 1) ~ "Scenario 1",
             (situation == "null" & targetdose == 3) ~ "Scenario 2",
             (situation == "null" & targetdose == 5) ~ "Scenario 3",
             (situation == "lowerDR" & targetdose == 1) ~ "Scenario 4",
             (situation == "lowerDR" & targetdose == 3) ~ "Scenario 5",
             (situation == "lowerDR" & targetdose == 5) ~ "Scenario 6",
             (situation == "higherDR" & targetdose == 1 & truedrprob_targetdose %in% c(65,55.1,65.1,65.2,55)) ~ "Scenario 7",
             (situation == "higherDR" & targetdose == 3) ~ "Scenario 9",
             TRUE ~ "Scenario 8"
           )) %>% 
  filter(scenario %in% c("Scenario 1","Scenario 2","Scenario 3","Scenario 5","Scenario 8","Scenario 9")) %>% 
  separate(measure, into=c("drdlt","component"),sep="_") %>% 
  mutate(dlttarget = as.character(dlttarget))
accuracy_pl$gen <- factor(accuracy_pl$gen,levels=c("Constant","Small Increasing","Small Decreasing","Large Decreasing"))

library(ggpubr)
library(grid)


############## AVERAGE OF ALL SCENARIOS

accuracy_pl_av <- accuracy_pl %>% group_by(dlttarget, n, gen, model,drdlt,component) %>% 
  summarize(valuemean = mean(value))

all_1 <- accuracy_pl_av %>% filter(component %in% c("bias","variance")) %>%  
  ggplot(aes(x=n,y=valuemean,group=interaction(dlttarget,model),color=model,linetype=dlttarget)) + 
  geom_line()+
  geom_point(size=0.8)+
  geom_hline(yintercept = 0,linetype=2,color="darkgrey")+
  facet_wrap(vars(gen,component,drdlt),nrow=4,scales = "free_y")+
  #coord_cartesian(ylim=c(0, 0.05))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",color='black',size=20,margin=margin(0,0,20,0)),
        #axis.title.x = element_text(hjust = 0.5,face="bold",color='black',size=15,margin=margin(20,0,0,0)),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA,fill=NA),
        legend.position="none")+
  scale_color_manual(values=c('#245D67','#7AC4D0'))+
  scale_x_continuous(name="", breaks=c(20,40,60,80,100))+
  labs(title="",xlab="Trial Sample Size (n)")


all <- annotate_figure(all_1,
                bottom = text_grob("Trial Sample Size", vjust=-0.2,hjust=0.5,size=13,face="bold"),
                right = text_grob(".       .",color="white"))+
  annotate("text", x = 0.26, y = 0.98, label = "Bias",size=5,fontface="bold")+
  annotate("text", x = 0.73, y = 0.98, label = "Variance",size=5,fontface="bold")+
  annotate("text", x = 0.16, y = 0.94, label = "DLT",size=3,color="#292929")+
  annotate("text", x = 0.37, y = 0.94, label = "DR",size=3,color="#292929")+
  annotate("text", x = 0.63, y = 0.94, label = "DLT",size=3,color="#292929")+
  annotate("text", x = 0.83, y = 0.94, label = "DR",size=3,color="#292929")+
  annotate("text", x=0.95,y=0.82,label="Constant",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.61,label="Small\nIncreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.41,label="Small\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.20,label="Large\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.15,y=0.03,label="40/20              50/30",size=3,color="#292929")+
  annotate("text", x=0.1,y=0.047,label="__       _ _",size=6,fontface="bold",color="#245D67")+
  annotate("text", x=0.82,y=0.03,label="Constant-Skeleton        Weibull-Skeleton",size=3,color="#292929")+
  annotate("text", x=0.626,y=0.061,label=".",size=20,fontface="bold",color="#245D67")+
  annotate("text", x=0.832,y=0.061,label=".",size=20,fontface="bold",color="#7AC4D0")
all
pdf("bv.pdf", width=6.5, height=7); all; dev.off()





mse_1 <- accuracy_pl_av %>% filter(component =="mse") %>%  
  ggplot(aes(x=n,y=valuemean,group=interaction(dlttarget,model),color=model,linetype=dlttarget)) + 
  geom_line()+
  geom_point(size=0.8)+
  geom_hline(yintercept = 0,linetype=2,color="darkgrey")+
  facet_wrap(vars(gen,component,drdlt),nrow=4,scales = "free_y")+
  #coord_cartesian(ylim=c(0, 0.05))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",color='black',size=20,margin=margin(0,0,20,0)),
        #axis.title.x = element_text(hjust = 0.5,face="bold",color='black',size=15,margin=margin(20,0,0,0)),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA,fill=NA),
        legend.position="none")+
  scale_color_manual(values=c('#245D67','#7AC4D0'))+
  scale_x_continuous(name="", breaks=c(20,40,60,80,100))+
  labs(title="",xlab="Trial Sample Size (n)")


mse <- annotate_figure(mse_1,
                       bottom = text_grob("Trial Sample Size", vjust=-0.2,hjust=0.45,size=13,face="bold"),
                       right = text_grob(".       .",color="white"))+
  annotate("text", x = 0.5, y = 0.98, label = "MSE",size=6,fontface="bold")+
  annotate("text", x = 0.27, y = 0.94, label = "DLT",size=4,color="#292929")+
  annotate("text", x = 0.73, y = 0.94, label = "DR",size=4,color="#292929")+
  annotate("text", x=0.95,y=0.82,label="Constant",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.61,label="Small\nIncreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.41,label="Small\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.20,label="Large\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.165,y=0.025,label="40/20              50/30",size=3,color="#292929")+
  annotate("text", x=0.1,y=0.042,label="__       _ _",size=6,fontface="bold",color="#245D67")+
  annotate("text", x=0.82,y=0.025,label="    Const-Skel       Weib-Skel",size=3,color="#292929")+
  annotate("text", x=0.67,y=0.052,label=".",size=20,fontface="bold",color="#245D67")+
  annotate("text", x=0.845,y=0.052,label=".",size=20,fontface="bold",color="#7AC4D0")
mse

pdf("mse.pdf", width=5, height=8); mse; dev.off()















###SCENARIO 2
s2_b <- accuracy_pl %>% filter(scenario=="Scenario 2",component=="bias") %>%  
  ggplot(aes(x=n,y=value,group=interaction(dlttarget,model),color=model,linetype=dlttarget)) + 
  geom_line()+
  geom_point(size=0.8)+
  geom_hline(yintercept = 0,linetype=2,color="darkgrey")+
  facet_wrap(vars(gen,component,drdlt),nrow=4,scales = "free_y")+
  #coord_cartesian(ylim=c(0, 0.05))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    #axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5,face="bold",color='black',size=20,margin=margin(0,0,20,0)),
    axis.title.x = element_text(hjust = 0.5,face="bold",color='black',size=15,margin=margin(20,0,0,0)),
    legend.background = element_blank(),
    legend.key = element_rect(color = NA,fill=NA),
    legend.position="none")+
  scale_color_manual(values=c('#245D67','#7AC4D0'))+
  scale_x_continuous(name="", breaks=c(20,40,60,80,100))+
  ggtitle("")


s2_v <- accuracy_pl %>% filter(scenario=="Scenario 2",component=="variance") %>%  
  ggplot(aes(x=n,y=value,group=interaction(dlttarget,model),color=model,linetype=dlttarget)) + 
  geom_line()+
  geom_point(size=0.8)+
  geom_hline(yintercept = 0,linetype=2,color="darkgrey")+
  facet_wrap(vars(gen,component,drdlt),nrow=4)+
  coord_cartesian(ylim=c(0, 0.03))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",color='black',size=20,margin=margin(0,0,20,0)),
        axis.title.x = element_text(hjust = 0.5,face="bold",color='black',size=15,margin=margin(20,0,0,0)),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA,fill=NA),
        legend.position="none")+
  scale_color_manual(values=c('#245D67','#7AC4D0'))+
  scale_x_continuous(name="", breaks=c(20,40,60,80,100))+
  scale_y_continuous(name="",breaks=c(0.00,0.02,0.04))+
  ggtitle("")

s2 <- ggarrange(s2_b ,s2_v,
          labels = NULL,
          ncol = 2,
          align = "hv", 
          font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

annotate_figure(s2, top = text_grob("Scenario 2", hjust = 0.4, size=20,face="bold"),
                bottom = text_grob("Trial Sample Size", vjust=-1.4,hjust=0.4,size=13,face="bold"),
                right = text_grob(".       .",color="white"))+
  annotate("text", x = 0.26, y = 0.92, label = "Bias",size=5,fontface="bold")+
  annotate("text", x = 0.73, y = 0.915, label = "Variance",size=5,fontface="bold")+
  annotate("text", x = 0.16, y = 0.88, label = "DLT",size=3,color="#292929")+
  annotate("text", x = 0.37, y = 0.88, label = "DR",size=3,color="#292929")+
  annotate("text", x = 0.63, y = 0.88, label = "DLT",size=3,color="#292929")+
  annotate("text", x = 0.83, y = 0.88, label = "DR",size=3,color="#292929")+
  annotate("text", x=0.96,y=0.76,label="Constant",size=3,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.59,label="Small\nIncreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.41,label="Small\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.24,label="Large\nDecreasing",size=3,color="#292929",angle="270")
  


##### SCENARIO 5
s5_b <- accuracy_pl %>% filter(scenario=="Scenario 5",component=="bias") %>%  
  ggplot(aes(x=n,y=value,group=interaction(dlttarget,model),color=model,linetype=dlttarget)) + 
  geom_line()+
  geom_point(size=0.8)+
  geom_hline(yintercept = 0,linetype=2,color="darkgrey")+
  facet_wrap(vars(gen,component,drdlt),nrow=4,scale="free_y")+
  #coord_cartesian(ylim=c(0, 0.07))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",color='black',size=20,margin=margin(0,0,20,0)),
        axis.title.x = element_text(hjust = 0.5,face="bold",color='black',size=15,margin=margin(20,0,0,0)),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA,fill=NA),
        legend.position="none")+
  scale_color_manual(values=c('#245D67','#7AC4D0'))+
  scale_x_continuous(name="", breaks=c(20,40,60,80,100))+
  ggtitle("")


s5_v <- accuracy_pl %>% filter(scenario=="Scenario 5",component=="variance") %>%  
  ggplot(aes(x=n,y=value,group=interaction(dlttarget,model),color=model,linetype=dlttarget)) + 
  geom_line()+
  geom_point(size=0.8)+
  geom_hline(yintercept = 0,linetype=2,color="darkgrey")+
  facet_wrap(vars(gen,component,drdlt),nrow=4)+
  coord_cartesian(ylim=c(0, 0.03))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",color='black',size=20,margin=margin(0,0,20,0)),
        axis.title.x = element_text(hjust = 0.5,face="bold",color='black',size=15,margin=margin(20,0,0,0)),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA,fill=NA),
        legend.position="none")+
  scale_color_manual(values=c('#245D67','#7AC4D0'))+
  scale_x_continuous(name="", breaks=c(20,40,60,80,100))+
  scale_y_continuous(name="",breaks=c(0.00,0.02,0.04))+
  ggtitle("")

s5 <- ggarrange(s5_b ,s5_v,
                labels = NULL,
                ncol = 2,
                align = "hv", 
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

annotate_figure(s5, top = text_grob("Scenario 5", hjust = 0.4, size=20,face="bold"),
                bottom = text_grob("Trial Sample Size", vjust=-1.4,hjust=0.4,size=13,face="bold"),
                right = text_grob(".       .",color="white"))+
  annotate("text", x = 0.26, y = 0.92, label ="Bias",size=5,fontface="bold")+
  annotate("text", x = 0.73, y = 0.915, label = "Variance",size=5,fontface="bold")+
  annotate("text", x = 0.16, y = 0.88, label = "DLT",size=3,color="#292929")+
  annotate("text", x = 0.37, y = 0.88, label = "DR",size=3,color="#292929")+
  annotate("text", x = 0.63, y = 0.88, label = "DLT",size=3,color="#292929")+
  annotate("text", x = 0.83, y = 0.88, label = "DR",size=3,color="#292929")+
  annotate("text", x=0.96,y=0.76,label="Constant",size=3,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.59,label="Small\nIncreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.41,label="Small\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.24,label="Large\nDecreasing",size=3,color="#292929",angle="270")



##### SCENARIO 9
s9_b <- accuracy_pl %>% filter(scenario=="Scenario 9",component=="bias") %>%  
  ggplot(aes(x=n,y=value,group=interaction(dlttarget,model),color=model,linetype=dlttarget)) + 
  geom_line()+
  geom_point(size=0.8)+
  geom_hline(yintercept = 0,linetype=2,color="darkgrey")+
  facet_wrap(vars(gen,component,drdlt),nrow=4,scale="free_y")+
  #coord_cartesian(ylim=c(0, 0.07))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",color='black',size=20,margin=margin(0,0,20,0)),
        axis.title.x = element_text(hjust = 0.5,face="bold",color='black',size=15,margin=margin(20,0,0,0)),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA,fill=NA),
        legend.position="none")+
  scale_color_manual(values=c('#245D67','#7AC4D0'))+
  scale_x_continuous(name="", breaks=c(20,40,60,80,100))+
  ggtitle("")


s9_v <- accuracy_pl %>% filter(scenario=="Scenario 9",component=="variance") %>%  
  ggplot(aes(x=n,y=value,group=interaction(dlttarget,model),color=model,linetype=dlttarget)) + 
  geom_line()+
  geom_point(size=0.8)+
  geom_hline(yintercept = 0,linetype=2,color="darkgrey")+
  facet_wrap(vars(gen,component,drdlt),nrow=4)+
  coord_cartesian(ylim=c(0, 0.03))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",color='black',size=20,margin=margin(0,0,20,0)),
        axis.title.x = element_text(hjust = 0.5,face="bold",color='black',size=15,margin=margin(20,0,0,0)),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA,fill=NA),
        legend.position="none")+
  scale_color_manual(values=c('#245D67','#7AC4D0'))+
  scale_x_continuous(name="", breaks=c(20,40,60,80,100))+
  scale_y_continuous(name="",breaks=c(0.00,0.02,0.04))+
  ggtitle("")

s9 <- ggarrange(s9_b ,s9_v,
                labels = NULL,
                ncol = 2,
                align = "hv", 
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

annotate_figure(s9, top = text_grob("Scenario 9", hjust = 0.4, size=20,face="bold"),
                bottom = text_grob("Trial Sample Size", vjust=-1.4,hjust=0.4,size=13,face="bold"),
                right = text_grob(".       .",color="white"))+
  annotate("text", x = 0.26, y = 0.92, label = "Bias",size=5,fontface="bold")+
  annotate("text", x = 0.73, y = 0.915, label = "Variance",size=5,fontface="bold")+
  annotate("text", x = 0.16, y = 0.88, label = "DLT",size=3,color="#292929")+
  annotate("text", x = 0.37, y = 0.88, label = "DR",size=3,color="#292929")+
  annotate("text", x = 0.63, y = 0.88, label = "DLT",size=3,color="#292929")+
  annotate("text", x = 0.83, y = 0.88, label = "DR",size=3,color="#292929")+
  annotate("text", x=0.96,y=0.76,label="Constant",size=3,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.59,label="Small\nIncreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.41,label="Small\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.24,label="Large\nDecreasing",size=3,color="#292929",angle="270")








###########  Bias, Variance, MSE together. No constant model. DLT and DR on same plot.
########### 


accuracy_pl_av <- accuracy_pl %>% group_by(dlttarget, n, gen, model,drdlt,component) %>% 
  summarize(valuemean = mean(value))
accuracy_pl_av$component <- factor(accuracy_pl_av$component,levels=c("bias","variance","mse"))

all_3 <- accuracy_pl_av %>% filter(model != "Constant-Skeleton") %>%  
  ggplot(aes(x=n,y=valuemean,group=interaction(dlttarget,drdlt),color=drdlt,linetype=dlttarget)) + 
  geom_line()+
  geom_point(size=0.8)+
  geom_hline(yintercept = 0,linetype=2,color="darkgrey")+
  facet_wrap(vars(gen,component),nrow=4,scales = "free_y")+
  #coord_cartesian(ylim=c(0, 0.05))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",color='black',size=20,margin=margin(0,0,5,0)),
        #axis.title.x = element_text(hjust = 0.5,face="bold",color='black',size=15,margin=margin(20,0,0,0)),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA,fill=NA),
        legend.position="none")+
  scale_color_manual(values=c('#245D67','#7AC4D0'))+
  scale_x_continuous(name="", breaks=c(20,40,60,80,100))+
  labs(title="",xlab="Trial Sample Size (n)")


all3 <- annotate_figure(all_3,
                       bottom = text_grob("Trial Sample Size", vjust=-0.2,hjust=0.5,size=13,face="bold"),
                       right = text_grob(".       .",color="white"))+
  annotate("text", x = 0.2, y = 0.98, label = "Bias",size=5,fontface="bold")+
  annotate("text", x = 0.49, y = 0.98, label = "Variance",size=5,fontface="bold")+
  annotate("text", x = 0.78, y = 0.98, label = "MSE",size=5,fontface="bold")+
  annotate("text", x=0.95,y=0.82,label="Constant",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.61,label="Small\nIncreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.41,label="Small\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.20,label="Large\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.15,y=0.03,label="40/20              50/30",size=3,color="#292929")+
  annotate("text", x=0.1,y=0.047,label="__       _ _",size=6,fontface="bold",color="#245D67")+
  annotate("text", x=0.82,y=0.03,label="DLT target         Dose Reduction target",size=3,color="#292929")+
  annotate("text", x=0.64,y=0.061,label=".",size=20,fontface="bold",color="#245D67")+
  annotate("text", x=0.778,y=0.061,label=".",size=20,fontface="bold",color="#7AC4D0")
all3
pdf("all3.pdf", width=6.5, height=7); all3; dev.off()




###########  Bias, Variance, MSE together. Only constant model. DLT and DR on same plot.
########### 


accuracy_pl_av <- accuracy_pl %>% group_by(dlttarget, n, gen, model,drdlt,component) %>% 
  summarize(valuemean = mean(value))
accuracy_pl_av$component <- factor(accuracy_pl_av$component,levels=c("bias","variance","mse"))

all_4 <- accuracy_pl_av %>% filter(model != "Weibull-Skeleton") %>%  
  ggplot(aes(x=n,y=valuemean,group=interaction(dlttarget,drdlt),color=drdlt,linetype=dlttarget)) + 
  geom_line()+
  geom_point(size=0.8)+
  geom_hline(yintercept = 0,linetype=2,color="darkgrey")+
  facet_wrap(vars(gen,component),nrow=4,scales = "free_y")+
  #coord_cartesian(ylim=c(0, 0.05))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",color='black',size=20,margin=margin(0,0,5,0)),
        #axis.title.x = element_text(hjust = 0.5,face="bold",color='black',size=15,margin=margin(20,0,0,0)),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA,fill=NA),
        legend.position="none")+
  scale_color_manual(values=c('#245D67','#7AC4D0'))+
  scale_x_continuous(name="", breaks=c(20,40,60,80,100))+
  labs(title="",xlab="Trial Sample Size (n)")


all4 <- annotate_figure(all_4,
                        bottom = text_grob("Trial Sample Size", vjust=-0.2,hjust=0.5,size=13,face="bold"),
                        right = text_grob(".       .",color="white"))+
  annotate("text", x = 0.2, y = 0.98, label = "Bias",size=5,fontface="bold")+
  annotate("text", x = 0.49, y = 0.98, label = "Variance",size=5,fontface="bold")+
  annotate("text", x = 0.78, y = 0.98, label = "MSE",size=5,fontface="bold")+
  annotate("text", x=0.95,y=0.82,label="Constant",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.61,label="Small\nIncreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.41,label="Small\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.95,y=0.20,label="Large\nDecreasing",size=3,color="#292929",angle="270")+
  annotate("text", x=0.15,y=0.03,label="40/20              50/30",size=3,color="#292929")+
  annotate("text", x=0.1,y=0.047,label="__       _ _",size=6,fontface="bold",color="#245D67")+
  annotate("text", x=0.82,y=0.03,label="DLT target         Dose Reduction target",size=3,color="#292929")+
  annotate("text", x=0.64,y=0.061,label=".",size=20,fontface="bold",color="#245D67")+
  annotate("text", x=0.778,y=0.061,label=".",size=20,fontface="bold",color="#7AC4D0")
all4
pdf("all4.pdf", width=6.5, height=7); all4; dev.off()



