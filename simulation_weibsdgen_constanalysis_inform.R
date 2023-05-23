library(rstan)
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

#read in Stan model analysis code. 
stancode_constant <- readLines("stancode_constant_inform.txt")
mymodel <- stan_model(model_code=stancode_constant)
#store event probability calculation functions
h_dr_const <- function(t,l1,l2,x0){
  (l1/(-l1-l2))*(exp((-l1*x0-l2*x0)*t)-1)
}
h_dltany_const <- function(t,l1,l2,b,x0,x1){
  (l1/(-l1-l2))*(exp((-l1*x0-l2*x0)*t)-1)-
    exp(-l2*(x1^exp(b))*t)*((l1*x0/(-l1*x0-l2*x0+l2*x1^exp(b)))*(exp((-l1*x0-l2*x0+l2*x1^exp(b))*t)-1))+
    (l2/(-l1-l2))*(exp((-l1*x0-l2*x0)*t)-1)
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
changegrid$type <- rep(c(rep("dec_small",6),rep("dec_small_DR35",6),rep("dec_small_DR65",6)),3)
changegrid$alpha <- 0.8
changegrid$lambda1 <- rep(c(0.0295,0.0214,0.0096,0.0206,0.0151,0.0068,0.0185,
                            0.0136,0.0061,0.0118,0.0086,0.0039,0.043,0.031,0.0142,0.032,0.023,0.0103),3)
changegrid$lambda2 <- rep(c(0.0145,0.01,0.0047,0.009,0.0063,0.0029,0.014,0.01,0.0046,0.0088,
                            0.0062,0.0028,0.015,0.01,0.0048,0.0095,0.0063,0.003),3)
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
    frame <- data.frame(matrix(nrow=sampsize, ncol=37+length(doses)*2))
    l1_percs <- paste("l1_",c("2.5","25","75","97.5"),sep='')
    l2_percs <- paste("l2_",c("2.5","25","75","97.5"),sep='')
    beta_percs <-paste("beta_",c("2.5","25","75","97.5"),sep='')
    dr0names <- paste("dr0prob_",as.character(doses[1:(length(doses)-1)]),sep="")
    dlt0names <- paste("dlt0prob_",as.character(doses[1:(length(doses)-1)]),sep="")
    names(frame) <- c("Patient","y1","delta1","y2","delta2","dose0","dose1","t1","t2","t3",
                      "accrualadd","accrualtime","y1_qum","y2_qum","current_delta1","current_delta2",
                      "current_y1","current_y2","l1","l2","beta",l1_percs,l2_percs,beta_percs, 
                      dr0names,dlt0names,"incstcy_indic","escalation","total_dlt","total_dr",
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
        Time1 <- frametemp$current_y1
        Time2 <- frametemp$current_y2
        Delt1 <- frametemp$current_delta1
        Delt2 <- frametemp$current_delta2
        Dose0 <- frametemp$dose0
        Dose1 <- frametemp$dose1
        #list data to go into stan model
        patientdata <- list(N=nrow(frametemp),y1=Time1,y2=Time2,delta1=Delt1,delta2=Delt2,dose0=Dose0,dose1=Dose1)
        #stan model
        fit <- sampling(mymodel,data=patientdata, 
                        iter=2000, chains=2,refresh=0)
        
        result <- rstan::extract(fit)
        #store model results
        l1 <- median(result$lambda1)
        l2 <- median(result$lambda2)
        b <- median(result$beta)
        frame[i,19:21]<- c(l1,l2,b)
        frame[i,22:25] <-  quantile(result$lambda1, c(0.025,0.25,0.75,0.975)) 
        frame[i,26:29] <-  quantile(result$lambda2, c(0.025,0.25,0.75,0.975)) 
        frame[i,30:33] <-  quantile(result$beta, c(0.025,0.25,0.75,0.975)) 
        #calculate and store probabilities
        drprobs <- h_dr_const(targettime,l1,l2,doses[1:(length(doses)-1)])
        dltprobs <- h_dltany_const(targettime,l1,l2,b,doses[1:(length(doses)-1)],doses[2:length(doses)])
        frame[i,34:(34+length(doses)-2)]<-drprobs
        frame[i,(34+length(doses)-1):(37+length(doses))]<-dltprobs
        #see if there are any monotonicity inconsistencies in calculated probabilites, DLT higher for smaller dose.
        if((frame[i,39]<frame[i,40])||
           (frame[i,40]<frame[i,41])||
           (frame[i,41]<frame[i,42])||
           (frame[i,42]<frame[i,43])){
          frame$incstcy_indic[i] <- 1
        }else{
          frame$incstcy_indic[i] <- 0
        }
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
        
        ############ Two different ways to choose next dose. 
        #choose closest dose to Dose Reduction Target and DLT target in total distance.
        dose0 <- doses[which((abs(drprobs-DRtarget)+abs(dltprobs-DLTtarget))==min((abs(drprobs-DRtarget)+abs(dltprobs-DLTtarget))))]
        dose1 <- doses[which(doses==dose0)+1]
        
        #choose dose set below DLT target and choose closest dose below DR target from that set. 
        #   doseDLTaccept <- doses[which(dltprobs<=DLTtarget)]
        #   doseDRaccept <- doses[which(drprobs<=DRtarget)]
        #   if(length(doseDLTaccept)>0 & length(doseDRaccept)>0){
        #     dose0=max(doseDRaccept[doseDRaccept %in% doseDLTaccept])
        #   }else(
        #     dose0 <- doses[order(doses)][2]
        #   )
        #   dose1 <- doses[which(doses==dose0)+1]
      }
    }
    store[[nsim]]<-frame
  }
  return(store)
}


library(parallel)
resultSet <- mclapply(c(1:54), runSim, N=3, mc.cores = numCores)
save(resultSet,file="/scratch/jeh6kw/results_weibsdgen_constanalysis_inform.Rdata")
Sys.time()-r1
