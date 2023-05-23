library(trialr)
library(codetools)
numCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))-1
options(mc.cores=numCores)
r1 <- Sys.time()


########################################################################################
#######################################################################################
########################################################################################

#Accrual Simulator

#function to generate a single independent observation based on starting dose with true model

gfunc_weib <- function(dose0,dose1,ralpha,rlambda1,rlambda2,rbeta0,rbeta1,gendoses){
  #one obs with DLT and DR times and indicators
  gendat <- data.frame(matrix(nrow=1,ncol=9))  
  names(gendat) <- c("y1","delta1","y2","delta2","dose0","dose1","t1","t2","t3")
  gendat$dose0 <- dose0
  gendat$dose1 <- dose1
  #how many weeks staying on trial total for an individual patient
  cens <- 52
  
  gendose0 <- gendoses[which(doses==dose0)]
  gendose1 <- gendoses[which(doses==dose1)]
  
  #generate the event times based on true model
  t1 <- (log(1-runif(1))/(-ralpha*rlambda1*exp(rbeta1*gendose0)))^(1/ralpha); gendat$t1 <- t1
  t2 <-(log(1-runif(1))/(-ralpha*rlambda2*exp(rbeta1*gendose0)))^(1/ralpha); gendat$t2 <- t2
  t3 <- (log(1-runif(1))/(-ralpha*rlambda1*exp(rbeta0+rbeta1*gendose1)))^(1/ralpha); gendat$t3 <- t3
  
  #record the event times with indicator notations. 
  if(t1>=cens & t2>=cens){
    gendat$y1<- cens
    gendat$y2 <- cens
    gendat$delta1 <- 0
    gendat$delta2 <- 0
  }
  
  if(t2<=t1 & t2<cens){
    gendat$delta1 <- 0
    gendat$delta2 <- 1
    gendat$y1 <- t2
    gendat$y2 <- t2
  }
  
  if(t2>t1 & t1<cens & (t3>=cens | t3<=t1)){
    gendat$delta1 <-1
    gendat$delta2 <-0
    gendat$y1 <- t1
    gendat$y2 <- cens
  }
  
  if(t2>t1 & t1<cens & t3<cens & t3>t1){
    gendat$delta1 <- 1
    gendat$delta2 <- 1
    gendat$y1 <-  t1
    gendat$y2 <- t3
  }
  return(gendat)
}


#simulation parameters
#static
accrualtime <- 0
accrualmean <- 10  #num of days average per patient accrual. Change to weeks in simulation.
doses <- c(1,0.65,0.4,0.2,0.1,0.05)
gendoses <- c(0.7,0.6,0.4,0.2,0.1,0.05)
dose0 <- doses[order(doses)][2]
dose1 <- doses[order(doses)][1]
targettime <- 52 #in weeks



#changing generator,type,sampsize,alpha,l1,l2,beta,beta0,beta1,DRtarget,DLTtarget
changegrid <- data.frame(matrix(nrow=54,ncol=10))
names(changegrid) <- c("generator","type","sampsize","alpha","lambda1","lambda2","beta0","beta1",
                       "DRtarget","DLTtarget")
changegrid$sampsize <- c(rep(101,18),rep(61,18),rep(31,18))
changegrid$generator <- "weibull"
changegrid$type <- rep(c(rep("inc_small",6),rep("dec_small",6),rep("dec_big",6)),3)
changegrid$alpha <- rep(c(rep(1.1,6),rep(0.8,6),rep(0.5,6)),3)
changegrid$lambda1 <- rep(c(0.012,0.0083,0.0038,0.0074,0.0053,0.0025,0.0295,
                            0.0214,0.0096,0.0206,0.0151,0.0068,0.08,0.058,0.026,0.06,0.044,0.0195),3)
changegrid$lambda2 <- rep(c(0.005,0.0033,0.0016,0.003,0.002,0.001,0.0145,0.01,0.0047,
                            0.009,0.0063,0.0029,0.043,0.03,0.014,0.028,0.02,0.009),3)
changegrid$beta0 <- -0.5
changegrid$beta1 <- 2
changegrid$DRtarget <- rep(c(rep(50,3),rep(40,3)),9)
changegrid$DLTtarget <- rep(c(rep(30,3),rep(20,3)),9)




####################
####################   Begin simulation

runSim <- function(n, N){
  print(n)
  generator <- changegrid$generator[n]
  type <- changegrid$type[n]
  sampsize <- changegrid$sampsize[n]
  ralpha <- changegrid$alpha[n]
  rlambda1 <- changegrid$lambda1[n]
  rlambda2 <- changegrid$lambda2[n]
  rbeta <- changegrid$beta[n]
  rbeta0 <- changegrid$beta0[n]
  rbeta1 <- changegrid$beta1[n]
  DRtarget <- changegrid$DRtarget[n]/100
  DLTtarget <- changegrid$DLTtarget[n]/100
  
  store <- vector(mode='list',length=N)
  
  
  for(nsim in 1:N){
    
    
    #set up results storage data frame for an individual simulation run. 
    frame <- data.frame(matrix(nrow=sampsize, ncol=24))
    names(frame) <- c("Patient","y1","delta1","y2","delta2","dose0","dose1","t1","t2","t3",
                      "accrualadd","accrualtime","y1_qum","y2_qum","current_delta1","current_delta2",
                      "current_y1","current_y2","incstcy_indic","escalation","total_dlt","total_dr",
                      "esc_after_DLT","deesc_after_non")
    
    
    #single simulation for loop
    for(i in 1:sampsize){
      #generate next patient observation. Either from constant generator or one of weibull generators
      temp <- gfunc_weib(dose0,dose1,ralpha,rlambda1,rlambda2,rbeta0,rbeta1,gendoses)
      
      
      
      #record accrual time and total time
      frame$Patient[i] <- i
      frame$accrualadd[i] <- rgeom(1,1/(accrualmean+1))/7
      frame$accrualtime[i] <- accrualtime
      accrualtime <- accrualtime+frame$accrualadd[i]
      
      #record parameters and event times
      frame$delta1[i] <- temp$delta1
      frame$delta2[i] <- temp$delta2
      frame$t1[i] <- temp$t1
      frame$t2[i] <- temp$t2
      frame$t3[i] <- temp$t3
      frame$y1[i] <- temp$y1
      frame$y2[i] <- temp$y2
      frame$dose0[i] <- temp$dose0
      frame$dose1[i] <- temp$dose1
      
      #record event times relative to beginning of trial
      frame$y1_qum[i] <- frame$y1[i] + frame$accrualtime[i]
      frame$y2_qum[i] <- frame$y2[i] + frame$accrualtime[i]
      
      #readjust censoring indicators for the patient's current state in trial, not final state
      for(j in 1:i){
        
        #truly censored obs already matured
        if(frame$y2_qum[j]-frame$accrualtime[i]<=0 & frame$delta1[j]==0 & frame$delta2[j]==0){
          frame$current_delta1[j] <- 0; frame$current_delta2[j] <- 0
          frame$current_y1[j] <- frame$y1[j]; frame$current_y2[j] <- frame$y2[j]
        }
        #truly DLT obs already matured
        if(frame$y2_qum[j]-frame$accrualtime[i]<=0 & frame$delta1[j]==0 & frame$delta2[j]==1){
          frame$current_delta1[j] <- 0; frame$current_delta2[j] <- 1
          frame$current_y1[j] <- frame$y1[j]; frame$current_y2[j] <- frame$y2[j]
        }
        #truly DR/cens obs already matured
        if(frame$y2_qum[j]-frame$accrualtime[i]<=0 & frame$delta1[j]==1 & frame$delta2[j]==0){
          frame$current_delta1[j] <- 1; frame$current_delta2[j] <- 0
          frame$current_y1[j] <- frame$y1[j]; frame$current_y2[j] <- frame$y2[j]
        }
        #truly DR/DLT obs already matured
        if(frame$y2_qum[j]-frame$accrualtime[i]<=0 & frame$delta1[j]==1 & frame$delta2[j]==1){
          frame$current_delta1[j] <- 1; frame$current_delta2[j] <- 1
          frame$current_y1[j] <- frame$y1[j]; frame$current_y2[j] <- frame$y2[j]
        }
        #truly any obs with no events matured yet
        if(frame$y1_qum[j]-frame$accrualtime[i]>0){
          frame$current_delta1[j] <- 0; frame$current_delta2[j] <- 0
          frame$current_y1[j] <- frame$accrualtime[i]-frame$accrualtime[j]+frame$accrualadd[i]; frame$current_y2[j] <- frame$accrualtime[i]-frame$accrualtime[j]+frame$accrualadd[i]
        }
        #truly DR/cens or DR/DLT with DR matured and cens/DLT not matured. 
        if(frame$y1_qum[j]-frame$accrualtime[i]<=0 & frame$y2_qum[j]-frame$accrualtime[i]>0){
          frame$current_delta1[j] <- 1; frame$current_delta2[j] <- 0
          frame$current_y1[j] <- frame$y1[j]; frame$current_y2[j] <- frame$accrualtime[i]+frame$accrualadd[i]-frame$y1_qum[j]+frame$y1[j]
        }
        
      }#end j loop
      
      #model doesn't like to fit with one obs, so 2 or more, fit model
      if(i>1){
        
        #make temporary data frame to plug into model
        frametemp <- frame[1:i,]
        
        ########
        ########   CODE TITECRM MODEL
        fit <- stan_crm(skeleton = c(0.05, 0.12, 0.20, 0.30, 0.50),
                        target = DLTtarget,
                        doses_given = match(frametemp$dose0,rev(doses))-1,
                        tox = frametemp$current_delta2,
                        weights = frametemp$current_y2/targettime,
                        model = 'empiric', beta_sd = sqrt(1.5), seed = 123, refresh=0)
        
        
        #Choose dose0 and dose1
        dose0 <- rev(doses)[fit$recommended_dose+1]
        dose1 <- rev(doses)[fit$recommended_dose]
        
        #escalation counter
        if(which(doses==frame$dose0[i])-which(doses==frame$dose0[i-1])==0){
          frame$escalation[i] <- 0}
        if(which(doses==frame$dose0[i])-which(doses==frame$dose0[i-1])==1){
          frame$escalation[i] <- -1}
        if(which(doses==frame$dose0[i])-which(doses==frame$dose0[i-1])==2){
          frame$escalation[i] <- -2}
        if(which(doses==frame$dose0[i])-which(doses==frame$dose0[i-1])==3){
          frame$escalation[i] <- -3}
        if(which(doses==frame$dose0[i])-which(doses==frame$dose0[i-1])==-1){
          frame$escalation[i] <- 1}
        if(which(doses==frame$dose0[i])-which(doses==frame$dose0[i-1])==-2){
          frame$escalation[i] <- 2}
        if(which(doses==frame$dose0[i])-which(doses==frame$dose0[i-1])==-3){
          frame$escalation[i] <- 3}
        #total tox and total reduction counter
        frame$total_dlt[i] <- sum(frame$current_delta2,na.rm = TRUE)
        frame$total_dr[i] <- sum(frame$current_delta1,na.rm = TRUE)
        
        #Check if dose escalates after period where a patient experienced a DLT, or deescalates after period with no DLTs
        frame$esc_after_DLT[i] <- ifelse((frame$total_dlt[i]-frame$total_dlt[i-1]>0)&(frame$escalation[i]>0),1,0)
        frame$deesc_after_non[i] <- ifelse((frame$total_dlt[i]-frame$total_dlt[i-1]==0)&(frame$escalation[i]<0),1,0)
        
        
        
      }
    }
    store[[nsim]]<-frame
  }
  return(store)
}


library(parallel)
resultSet <- mclapply(c(1:54), runSim, N=3, mc.cores = numCores)
save(resultSet,file="/scratch/jeh6kw/results_titecrm_weibullgen.Rdata")
Sys.time()-r1
