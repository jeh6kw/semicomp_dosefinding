#weib bd gen

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
  t1 <- (log(1-runif(1))/(-rlambda1*exp(rbeta1*gendose0)))^(1/ralpha); gendat$t1 <- t1
  t2 <- (log(1-runif(1))/(-rlambda2*exp(rbeta1*gendose0)))^(1/ralpha); gendat$t2 <- t2
  t3 <- (log(1-runif(1))/(-rlambda2*exp(rbeta0+rbeta1*gendose1)))^(1/ralpha); gendat$t3 <- t3
  
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
sm <- stan_model(file="empiric_stan2.stan", model_name='prob2',verbose=FALSE)
titecrmmc <- function(x,doselabel,y,follow,alpha0,Tmax,target1,target2){
  # x dose level
  # doselabel    numbers used in regression
  # y grade
  # weight
  # alpha0 intercept
  pos_util <- which(follow!=0.1)
  
  if (length(pos_util)>0){
    N=length(pos_util)
    K=length(doselabel)
    #x1p <- doselabel[x]
    w1p <- w2p <- follow/Tmax
    w1p[y==2] <- rep(1, sum(y==2))
    w2p[y==3] <- rep(1, sum(y==3))
    
    data_p = list( a=alpha0, N=N, K=K, x=array(x[pos_util]), y=array(y[pos_util]), 
                   w1=array(w1p[pos_util]), w2=array(w2p[pos_util]), d=doselabel)
    
    fit <- sampling(sm, data=data_p ,iter=4000, cores=1, chains=4, refresh=0, control = list(adapt_delta = 0.8))
    #coeff=get_posterior_mean(fit)
    
    est <- extract(fit, pars="theta")
    theta <- apply(est$theta,2,median)
    theta <- min(theta)
    
    est2 <- extract(fit, pars=c("beta", "gam"))
    bet <- exp(median(est2$beta))
    gam <- exp(median(est2$gam))
    
    #p1tox<-pnorm(alpha0+est[1]*doselabel)
    #p2tox<-pnorm(alpha0+est[1]*doselabel-est[2])
    #cur1<-which(abs(est[3:(2+K),5]-target1)==min(abs(est[3:(2+K),5]-target1)))
    #cur2<-which(abs(est[(3+K):(2+2*K),5]-target2)==min(abs(est[(3+K):(2+2*K),5]-target2)))
    cur<-order(abs(doselabel-theta))[1]
  } else {
    cur=max(x)
    bet=-1
    gam=-1
  }
  
  
  list(newdose=cur, p1tox=doselabel^bet, p2tox=doselabel^(gam+bet))
  
}


#simulation parameters
#static

accrualmean <- 10  #num of days average per patient accrual. Change to weeks in simulation.
doses <- c(1,0.65,0.4,0.2,0.1,0.05)
gendoses <- c(0.7,0.6,0.4,0.2,0.1,0.05)
dose0 <- doses[order(doses)][2]
dose1 <- doses[order(doses)][1]
targettime <- 52 #in weeks
J <- 12
doselabel <- c(0.14005, 0.25000, 0.37620, 0.50185, 0.61495)
alpha0 = 1
K=5


#changing generator,type,sampsize,alpha,l1,l2,beta,beta0,beta1,DRtarget,DLTtarget
changegrid <- data.frame(matrix(nrow=54,ncol=10))
names(changegrid) <- c("generator","type","sampsize","alpha","lambda1","lambda2","beta0","beta1",
                       "DRtarget","DLTtarget")
changegrid$sampsize <- c(rep(101,18),rep(61,18),rep(31,18))
changegrid$generator <- "weibull"
changegrid$type <- rep(c(rep("decbig",6),rep("decbig_DR35",6),rep("decbig_DR65",6)),3)
changegrid$alpha <- 0.5
changegrid$lambda1 <- rep(c(0.104,0.057,0.031,0.068,0.037,0.0202,0.0613,0.034,0.0185,0.0374,0.0204,0.0112,0.167,0.093,0.0504,0.1085,0.06,0.0326),3)
changegrid$lambda2 <- rep(c(0.049,0.028,0.015,0.029,0.0165,0.0087,0.0454,0.026,0.0139,0.0275,0.0153,0.0082,0.052,0.031,0.0161,0.031,0.0178,0.0093),3)
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
    accrualtime <- 0
    x <- c(1,1)
    
    #single simulation for loop
    for(i in 1:sampsize){
      
      #generate next patient observation. Either from constant generator or one of weibull generators
      temp <- gfunc_weib(dose0,dose1,ralpha,rlambda1,rlambda2,rbeta0,rbeta1,gendoses)
      
      
      
      #record accrual time and total time
      frame$Patient[i] <- i
      frame$accrualadd[i] <- max(rgeom(1,1/(accrualmean+1))/7,0.1)
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
        
        eventobs <- vector(length=nrow(frametemp))
        cycleobs <- vector(length=nrow(frametemp))
        
        
        
        for(ob in 1:nrow(frametemp)){
          if(frametemp$current_delta1[ob]==0 & frametemp$current_delta2[ob]==0){eventobs[ob]<-1}
          if(frametemp$current_delta1[ob]==1 & frametemp$current_delta2[ob]==0){eventobs[ob]<-2}
          if(frametemp$current_delta2[ob]==1){eventobs[ob]<-3}
          
        }
        
        cycleobs <- ifelse((frametemp$accrualtime[nrow(frametemp)]-frametemp$accrualtime+frametemp$accrualadd[nrow(frametemp)])/(52/12)>12,12,(frametemp$accrualtime[nrow(frametemp)]-frametemp$accrualtime+frametemp$accrualadd[nrow(frametemp)])/(52/12))
        
        data_reg <- cbind(frametemp$Patient,cycleobs,eventobs)
        
        
        results <- suppressWarnings(titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,DRtarget,DLTtarget))
        closeAllConnections()
        newdose <- min(results$newdose, max(x)+1, K)
        
        dose0 <- rev(doses)[newdose+1]
        dose1 <- rev(doses)[newdose]
        x <- c(x,newdose)
        
      }
    }
    
    eventobs_full <- vector(length=nrow(frametemp))
    cycleobs_full <- vector(length=nrow(frametemp))
    
    for(ob in 1:nrow(frametemp)){
      if(frame$delta1[ob]==0 & frame$delta2[ob]==0){eventobs_full[ob]<-1}
      if(frame$delta1[ob]==1 & frame$delta2[ob]==0){eventobs_full[ob]<-2}
      if(frame$delta2[ob]==1){eventobs_full[ob]<-3}
    }
    
    cycleobs_full <- rep(J,nrow(frame))
    
    data_reg_full <- cbind(frame$Patient,cycleobs_full,eventobs_full)
    
    results_full <- suppressWarnings(titecrmmc(x,doselabel,y=data_reg_full[,3],follow=data_reg_full[,2],alpha0,Tmax=J,DRtarget,DLTtarget))
    closeAllConnections()
    newdose_full <- min(results$newdose, max(x)+1, K)
    
    dose0 <- rev(doses)[newdose_full+1]
    
    frame$CRMMCfullinfodose[nrow(frame)] <- dose0
    
    
    store[[nsim]]<-frame
  }
  return(store)
}



library(parallel)
resultSet <- mclapply(c(1:54), runSim, N=100, mc.cores = numCores)
save(resultSet,file="/home/jeh6kw/results_weibbdgen_titecrm_mcext.Rdata")
Sys.time()-r1