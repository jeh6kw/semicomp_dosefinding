########Patients followed max of 52 weeks no matter what. 
library(trialr)
load("results_titecrm_constantgen.Rdata")

store1 <- data.frame(matrix(nrow=18,ncol=13))
names(store1) <- c("generator","type","sampsize","lambda1","lambda2","beta","DRtarget","DLTtarget","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec")
store1$sampsize <- c(rep(101,6),rep(61,6),rep(31,6))
store1$generator <- "constant"
store1$type <- "null"
store1$lambda1 <- rep(c(0.162,0.042,0.017,0.111,0.029,0.0113),3)
store1$lambda2 <- rep(c(0.058,0.0175,0.0074,0.038,0.011,0.0045),3)
store1$beta <- -0.5
store1$DRtarget <- rep(c(rep(50,3),rep(40,3)),3)
store1$DLTtarget <- rep(c(rep(30,3),rep(20,3)),3)
doses <- c(1,0.65,0.4,0.2,0.1,0.05)

#### CRM after 52 - Constant scenarios

#Calculate cut off dose from the last estimates on each simulated trial
for(i in 1:nrow(store1)){
  print(i)
  #counter to pick out the correct row of the results_titecrm_constantgen.Rdata
  DLTtarget <- store1$DLTtarget[i]/100
  DRtarget <- store1$DRtarget[i]/100
  s1 <- 0;  s2 <- 0;  s3 <- 0;  s4 <- 0;  s5 <- 0
  
  
    
     
  sett <- resultSet[[i]]
    
  for(m in 1:1000){
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
  store1[i,9:13] <- c(s1,s2,s3,s4,s5)/10
  
}
write.csv(store1,"crmafter52_const.csv")





