
#### CRM after 52 - Weibull scenarios
library(trialr)
#changing generator,type,sampsize,alpha,l1,l2,beta,beta0,beta1,DRtarget,DLTtarget
store2 <- data.frame(matrix(nrow=54,ncol=15))
names(store2) <- c("generator","type","sampsize","alpha","lambda1","lambda2","beta0","beta1",
                   "DRtarget","DLTtarget","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec")
store2$sampsize <- c(rep(101,18),rep(61,18),rep(31,18))
store2$generator <- "weibull"
store2$type <- rep(c(rep("inc_small",6),rep("dec_small",6),rep("dec_big",6)),3)
store2$alpha <- rep(c(rep(1.1,6),rep(0.8,6),rep(0.5,6)),3)
store2$lambda1 <- rep(c(0.012,0.0083,0.0038,0.0074,0.0053,0.0025,0.0295,
                        0.0214,0.0096,0.0206,0.0151,0.0068,0.08,0.058,0.026,0.06,0.044,0.0195),3)
store2$lambda2 <- rep(c(0.005,0.0033,0.0016,0.003,0.002,0.001,0.0145,0.01,0.0047,
                        0.009,0.0063,0.0029,0.043,0.03,0.014,0.028,0.02,0.009),3)
store2$beta0 <- -0.5
store2$beta1 <- 2
store2$DRtarget <- rep(c(rep(50,3),rep(40,3)),9)
store2$DLTtarget <- rep(c(rep(30,3),rep(20,3)),9)
doses <- c(1,0.65,0.4,0.2,0.1,0.05)
############################################################################################################################
############################################################################################################################



#Calculate cut off dose from the last estimates on each simulated trial
for(i in 1:nrow(store2)){
  print(i)
  #counter to pick out the correct row of the results_titecrm_constantgen.Rdata
  DLTtarget <- store2$DLTtarget[i]/100
  DRtarget <- store2$DRtarget[i]/100
  s1 <- 0;  s2 <- 0;  s3 <- 0;  s4 <- 0;  s5 <- 0
  
  
  
  for(j in 1:10){
    #pull in right data set
    load(paste("results_titecrm_weibullgen_",j,".Rdata",sep=''))
    sett <- resultSet[[i]]
    
    for(m in 1:100){
      crnt <- sett[[m]]
      fit <- stan_crm(skeleton = c(0.05, 0.12, 0.20, 0.30, 0.50),
                      target = DLTtarget,
                      doses_given = match(crnt$dose0[1:(nrow(crnt)-1)],rev(doses))-1,
                      tox = crnt$delta2[1:(nrow(crnt)-1)],
                      model = 'empiric', beta_sd = sqrt(1.5), seed = 123, refresh=0)
      
      
      #Choose dose0 and dose1
      dose0 <- rev(doses)[fit$recommended_dose+1]
      dose1 <- rev(doses)[fit$recommended_dose]
      
      if(dose0==0.1){s1 <- s1+1}
      if(dose0==0.2){s2 <- s2+1}
      if(dose0==0.4){s3 <- s3+1}
      if(dose0==0.65){s4 <- s4+1}
      if(dose0==1){s5 <- s5+1}
    }#end m
    rm(list=paste("resultSet"))
  } #end j
  store2[i,11:15] <- c(s1,s2,s3,s4,s5)/10
  
}
write.csv(store2,"crmafter52_weibull.csv")