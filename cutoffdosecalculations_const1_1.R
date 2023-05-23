#Constant Analysis

#Build results table. Include DLTtarget, DRtarget, generator, trueDLT dose, DR scenario
store <- data.frame(matrix(nrow=54*8,ncol=11))
names(store) <- c("dlttarget","drtarget","targetdose","n","DRscenario","generator","dose1_selec","dose2_selec","dose3_selec","dose4_selec","dose5_selec")
store$dlttarget <- rep(rep(c(rep(30,3),rep(20,3)),9),8)
store$drtarget <- rep(rep(c(rep(50,3),rep(40,3)),9),8)
store$targetdose <- rep(rep(c(1,3,5),18),8)
store$n <- rep(c(rep(100,18),rep(60,18),rep(30,18)),8)
store$DRscenario <- rep(c(rep("null",6),rep("lowerDR",6),rep("higherDR",6)),8)
store$generator <- c(rep("const_uninf",54),rep("const_inf",54),rep("weibsd_uninf",54),rep("weibsd_inf",54),
                     rep("weibinc_inf",54),rep("weibinc_uninf",54),rep("weibbd_inf",54),rep("weibbd_uninf",54))
doses <- c(1,0.65,0.4,0.2,0.1,0.05)

############################################################################################################################
############################################################################################################################



cup <- 0
#Calculate cut off dose from the last estimates on each simulated trial
for(i in 55:nrow(store)){
  print(i)
  #counter to pick out the correct row of each of the 10 identical .Rdata's for each situation
  cup <- cup + 1
  if(cup>54){cup <- 1}
  DLTtarget <- store$dlttarget[i]/100
  DRtarget <- store$drtarget[i]/100
  s1 <- 0;  s2 <- 0;  s3 <- 0;  s4 <- 0;  s5 <- 0
  
  #set up correct data set extension
  if(i<=54){ext <- "results_constgen_constanalysis_uninform_"   }
  if(54<i && i<=108){ext <- "results_constgen_constanalysis_inform_" }
  if(108<i && i<=162){ext <- "results_weibsdgen_constanalysis_uninform_"  }
  if(162<i && i<=216){ext <-  "results_weibsdgen_constanalysis_inform_" }
  if(216<i && i<=270){ext <-  "results_weibincgen_constanalysis_inform_" }
  if(270<i && i<=324){ext <-  "results_weibincgen_constanalysis_uninform_" }
  if(324<i && i<=378){ext <-  "results_weibbdgen_constanalysis_inform_" }
  if(378<i && i<=432){ext <-  "results_weibbdgen_constanalysis_uninform_" }
  
  
  
  for(j in 1:10){
    #pull in right data set
    load(paste(ext,j,".Rdata",sep=''))
    resultSet <- get(paste("resultSet",j,sep='')) 
    sett <- resultSet[[cup]]
    
    for(m in 1:100){
      crnt <- sett[[m]]
      lastrow <- crnt[nrow(crnt)-1,]
      dltprobs <- lastrow[startsWith(names(lastrow),"dlt0")]
      drprobs <- lastrow[startsWith(names(lastrow),"dr0")]
      #choose dose set below DLT target and choose closest dose below DR target from that set.
      doseDLTaccept <- doses[which(dltprobs<=DLTtarget)]
      doseDRaccept <- doses[which(drprobs<=DRtarget)]
      if(length(doseDLTaccept)>0 & length(doseDRaccept)>0){  
        dose0=max(doseDRaccept[doseDRaccept %in% doseDLTaccept])
      }else{dose0 <- doses[order(doses)][2]} 
      dose1 <- doses[which(doses==dose0)+1]
      
      if(dose0==0.1){s1 <- s1+1}
      if(dose0==0.2){s2 <- s2+1}
      if(dose0==0.4){s3 <- s3+1}
      if(dose0==0.65){s4 <- s4+1}
      if(dose0==1){s5 <- s5+1}
      
      
      
    }#end m
    rm(list=paste("resultSet",j,sep=''))
  }#end j
  store$dose1_selec[i] <- s1/10
  store$dose2_selec[i] <- s2/10
  store$dose3_selec[i] <- s3/10
  store$dose4_selec[i] <- s4/10
  store$dose5_selec[i] <- s5/10
  
}
write.csv(store,"cutoffselection.csv")









