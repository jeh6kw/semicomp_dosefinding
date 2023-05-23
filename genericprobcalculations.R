library(tidyverse)
library(ggrepel)
library(cowplot)
library(ggtext)
#################
## Create generic, numerical integral function to plug in any hazard functions
##
rdose <- c(0.05,0.1,0.2,0.4,0.65)
idose <-  c(0.1,0.2,0.4,0.65,1)

rdose <- c(0.05,0.1,0.2,0.4,0.6)
idose <- c(0.1,0.2,0.4,0.6,0.7)

# #Change h's as necessary, integrand functions don't change
# #Skeleton-Constant 
# 
h1_const <- function(r,x0){  lambda1*x0  }
h2_const <- function(r,x0){  lambda2*x0  }
h3_const <- function(r,x1){  lambda2*x1^exp(beta)  }


#Cox-Weibull
#decreasing - extreme
#decreasing - moderate
#increasing - moderate
h1_weib <- function(r,x0){  alpha*lambda1*r^(alpha-1)*exp(beta1*x0)  }
h2_weib <- function(r,x0){  alpha*lambda2*r^(alpha-1)*exp(beta1*x0)  }
h3_weib <- function(r,x1){  alpha*lambda2*r^(alpha-1)*exp(beta0+beta1*x1)  }


h_12_integrand_const <- function(r,x0){
  h1_const(r,x0)*exp(-lambda1*x0*r)*exp(-lambda2*x0*r)
}
h_13_integrand_const <- function(r,x0){
  h2_const(r,x0)*exp(-lambda1*x0*r)*exp(-lambda2*x0*r)
}
h_123_integrand_const <- function(r,x0,x1,t){
  h1_const(r,x0)*exp(-lambda1*x0*r)*exp(-lambda2*x0*r)*(1-exp(-lambda2*t*x1^exp(beta))/exp(-lambda2*r*x1^exp(beta)))
}

h_12_integrand_weib <- function(r,x0){
  h1_weib(r,x0)*exp(-lambda1*r^alpha*exp(beta1*x0))*exp(-lambda2*r^alpha*exp(beta1*x0))
}
h_13_integrand_weib <- function(r,x0){
  h2_weib(r,x0)*exp(-lambda1*r^alpha*exp(beta1*x0))*exp(-lambda2*r^alpha*exp(beta1*x0))
}
h_123_integrand_weib <- function(r,x0,x1,t){
  h1_weib(r,x0)*exp(-lambda1*r^alpha*exp(beta1*x0))*exp(-lambda2*r^alpha*exp(beta1*x0))*(1-exp(-lambda2*t^alpha*exp(beta0+beta1*x1))/exp(-lambda2*r^alpha*exp(beta0+beta1*x1)))
}

#Endpoint Probs
targettime=52

p_12 <- function(targettime,x0){integrate(h_12_integrand_weib,lower=0,upper=targettime,x0=x0)[[1]]} 
p_13 <- function(targettime,x0){integrate(h_13_integrand_weib,lower=0,upper=targettime,x0=x0)[[1]]}
#need to keep it as targettime. Used inside function 
p_13_123 <- function(targettime,x0,x1){integrate(h_123_integrand_weib,lower=0,upper=targettime,x0=x0,x1=x1,t=targettime)[[1]]+integrate(h_13_integrand_weib,lower=0,upper=targettime,x0=x0)[[1]]}

alpha= 0.5;lambda1= 0.0326;lambda2= 0.0093; beta0 = -0.5; beta1 = 2; 

for(i in 1:length(rdose)){
  x0 <- idose[i]
  x1 <- rdose[i]
  a <- p_12(targettime,x0)
  b <- p_13(targettime,x0)
  c <- p_13_123(targettime,x0,x1)
  print(c(a,b,c))
}


#### Plot Hazards
#Constant
h1 <- function(r,x0){  lambda1*x0  }
h2 <- function(r,x0){  lambda2*x0  }
h3 <- function(r,x1){  lambda2*x1^exp(beta)  }
lambda1 = 0.042; lambda2 = 0.0175; beta = -1;
idose <-  c(0.1,0.2,0.4,0.65,1)
rdose <- c(0.05,0.1,0.2,0.4,0.65)
cdat <- data.frame(matrix(nrow=52*5,ncol=9)); names(cdat) <- c("lambda1","lambda2","beta","time","idose","rdose","h1","h2","h3")
cdat$time <- rep(seq(from=1, to=52,by=1),5)
cdat$lambda1 <- 0.042; cdat$lambda2 <- 0.0175; cdat$beta <- -0.5; cdat$idose <- rep(idose,52); cdat$rdose <- rep(rdose,52)

for(i in 1:52){
  for(j in 1:5){
    x1 <- rdose[j]
    x0 <- idose[j]
    cdat$h3[(i-1)*5+j] <- h3(i,x1)
    cdat$h2[(i-1)*5+j] <- h2(i,x0)
    cdat$h1[(i-1)*5+j] <- h1(i,x0)
  }
}
cdat_long <- pivot_longer(cdat,cols=starts_with("h"),names_to = "DRorDLT2",values_to = "HazardRate")
cdat_long$doseandhaz <- paste(cdat_long$idose,cdat_long$DRorDLT2)
cdat_long$label <- NA
cdat_long$label[which(cdat_long$time == max(cdat_long$time ))] <- cdat_long$DRorDLT2[which(cdat_long$time  == max(cdat_long$time ))]


w1 <- cdat_long %>% filter(idose==0.4) %>%  
  ggplot(aes(x=time,y=HazardRate,color=DRorDLT2,group=doseandhaz))+
  geom_line(size=2)+
  geom_label_repel(aes(label = label),nudge_x=1,na.rm=TRUE)+
  scale_color_manual(values=c("#9cb99d","#245d67","#003046"))+
  coord_cartesian(ylim = c(0, 0.055))+
  ylab("Hazard Rate")+
  xlab("Weeks")+
  theme(panel.grid.major.y =element_line(colour="#404040",linetype=2),
        panel.background = element_rect(fill = "white"),
        title= element_text("constant"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  annotate(geom="text", x=45, y=0.05, label="Constant",
           color="#404040",size=6)




#Weibull increasing
h1 <- function(r,x0){  alpha*lambda1*r^(alpha-1)*exp(beta1*x0)  }
h2 <- function(r,x0){  alpha*lambda2*r^(alpha-1)*exp(beta1*x0)  }
h3 <- function(r,x1){  alpha*lambda2*r^(alpha-1)*exp(beta0+beta1*x1)  }
alpha = 1.1; lambda1 = 0.006; lambda2 = 0.0028; beta0 = -0.5; beta1 = 2
idose <- c(0.1,0.2,0.4,0.6,0.7)
rdose <- c(0.05,0.1,0.2,0.4,0.6)

wincdat <- data.frame(matrix(nrow=52*5,ncol=11)); names(wincdat) <- c("alpha","lambda1","lambda2","beta0","beta1","time","idose","rdose","h1","h2","h3")
wincdat$time <- rep(seq(from=1, to=52,by=1),5)
wincdat$alpha <- 1.1; wincdat$lambda1 <- 0.006; wincdat$lambda2 <-0.0028; wincdat$beta0 <- -0.5; wincdat$beta1 <- 2; 
wincdat$idose <- c(rep(idose[1],52),rep(idose[2],52),rep(idose[3],52),rep(idose[4],52),rep(idose[5],52))
wincdat$rdose <- c(rep(rdose[1],52),rep(rdose[2],52),rep(rdose[3],52),rep(rdose[4],52),rep(rdose[5],52))

for(i in 1:nrow(wincdat)){
  x0 <- wincdat$idose[i]
  x1 <- wincdat$rdose[i]
  wincdat$h1[i] <- h1(wincdat$time[i],x0)
  wincdat$h2[i] <- h2(wincdat$time[i],x0)
  wincdat$h3[i] <- h3(wincdat$time[i],x1)
}
  
wincdat_long <- pivot_longer(wincdat,cols=starts_with("h"),names_to = "DRorDLT2",values_to = "HazardRate")
wincdat_long$doseandhaz <- paste(wincdat_long$idose,wincdat_long$DRorDLT2)
wincdat_long$label <- NA
wincdat_long$label[which(wincdat_long$time == max(wincdat_long$time ))] <- wincdat_long$DRorDLT2[which(wincdat_long$time  == max(wincdat_long$time ))]


w2 <- wincdat_long %>% filter(idose==0.4) %>% 
  ggplot(aes(x=time,y=HazardRate,color=DRorDLT2,group=doseandhaz))+
  geom_line(size=2)+
  geom_label_repel(aes(label = label),nudge_x=1,na.rm=TRUE)+
  scale_color_manual(values=c("#9cb99d","#245d67","#003046"))+
  coord_cartesian(ylim = c(0, 0.055))+
  ylab("Hazard Rate")+
  xlab("Weeks")+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  annotate(geom="text", x=40, y=0.05, label="Small Increasing",
           color="#404040",size=6)




#Weibull decreasing small
h1 <- function(r,x0){  alpha*lambda1*r^(alpha-1)*exp(beta1*x0)  }
h2 <- function(r,x0){  alpha*lambda2*r^(alpha-1)*exp(beta1*x0)  }
h3 <- function(r,x1){  alpha*lambda2*r^(alpha-1)*exp(beta0+beta1*x1)  }
alpha = .8; lambda1 = 0.0146; lambda2 = 0.0077; beta0 = -0.5; beta1 = 2
idose <- c(0.1,0.2,0.4,0.6,0.7)
rdose <- c(0.05,0.1,0.2,0.4,0.6)
wdecsdat <- data.frame(matrix(nrow=52*5,ncol=11)); names(wdecsdat) <- c("alpha","lambda1","lambda2","beta0","beta1","time","idose","rdose","h1","h2","h3")
wdecsdat$time <- rep(seq(from=1, to=52,by=1),5)
wdecsdat$alpha <- .8; wdecsdat$lambda1 <- 0.0146; wdecsdat$lambda2 <-0.0077; wdecsdat$beta0 <- -0.5; wdecsdat$beta1 <- 2; 
wdecsdat$idose <- c(rep(idose[1],52),rep(idose[2],52),rep(idose[3],52),rep(idose[4],52),rep(idose[5],52))
wdecsdat$rdose <- c(rep(rdose[1],52),rep(rdose[2],52),rep(rdose[3],52),rep(rdose[4],52),rep(rdose[5],52))

for(i in 1:nrow(wdecsdat)){
  x0 <- wdecsdat$idose[i]
  x1 <- wdecsdat$rdose[i]
  wdecsdat$h1[i] <- h1(wdecsdat$time[i],x0)
  wdecsdat$h2[i] <- h2(wdecsdat$time[i],x0)
  wdecsdat$h3[i] <- h3(wdecsdat$time[i],x1)
}
wdecsdat_long <- pivot_longer(wdecsdat,cols=starts_with("h"),names_to = "DRorDLT2",values_to = "HazardRate")
wdecsdat_long$doseandhaz <- paste( wdecsdat_long$idose, wdecsdat_long$DRorDLT2)
wdecsdat_long$label <- NA
wdecsdat_long$label[which(wdecsdat_long$time == max(wdecsdat_long$time ))] <- wdecsdat_long$DRorDLT2[which(wdecsdat_long$time  == max(wdecsdat_long$time ))]


w3 <- wdecsdat_long %>% filter(idose==0.4) %>% 
  ggplot(aes(x=time,y=HazardRate,color=DRorDLT2,group=doseandhaz))+
  geom_line(size=2)+
  geom_label_repel(aes(label = label),nudge_x=1,na.rm=TRUE)+
  coord_cartesian(ylim = c(0, 0.055))+
  scale_color_manual(values=c("#9cb99d","#245d67","#003046"))+
  ylab("Hazard Rate")+
  xlab("Weeks")+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),
        panel.background = element_rect(fill = "white"),
        title= element_text("increasing"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
    annotate(geom="text", x=40, y=0.05, label="Small Decreasing",
         color="#404040",size=6)


#Weibull decreasing large
h1 <- function(r,x0){  alpha*lambda1*r^(alpha-1)*exp(beta1*x0)  }
h2 <- function(r,x0){  alpha*lambda2*r^(alpha-1)*exp(beta1*x0)  }
h3 <- function(r,x1){  alpha*lambda2*r^(alpha-1)*exp(beta0+beta1*x1)  }
alpha = .5; lambda1 = 0.039; lambda2 = 0.022; beta0 = -0.5; beta1 = 2
idose <- c(0.1,0.2,0.4,0.6,0.7)
rdose <- c(0.05,0.1,0.2,0.4,0.6)


wdeclardat <- data.frame(matrix(nrow=52*5,ncol=11)); names(wdeclardat) <- c("alpha","lambda1","lambda2","beta0","beta1","time","idose","rdose","h1","h2","h3")
wdeclardat$time <- rep(seq(from=1, to=52,by=1),5)
wdeclardat$alpha <- .5; wdeclardat$lambda1 <- 0.039; wdeclardat$lambda2 <-0.022; wdeclardat$beta0 <- -0.5; wdeclardat$beta1 <- 2; 
wdeclardat$idose <- c(rep(idose[1],52),rep(idose[2],52),rep(idose[3],52),rep(idose[4],52),rep(idose[5],52))
wdeclardat$rdose <- c(rep(rdose[1],52),rep(rdose[2],52),rep(rdose[3],52),rep(rdose[4],52),rep(rdose[5],52))

for(i in 1:nrow(wdeclardat)){
  x0 <- wdeclardat$idose[i]
  x1 <- wdeclardat$rdose[i]
  wdeclardat$h1[i] <- h1(wdeclardat$time[i],x0)
  wdeclardat$h2[i] <- h2(wdeclardat$time[i],x0)
  wdeclardat$h3[i] <- h3(wdeclardat$time[i],x1)
}
wdeclardat_long <- pivot_longer(wdeclardat,cols=starts_with("h"),names_to = "DRorDLT2",values_to = "HazardRate")
wdeclardat_long$doseandhaz <- paste( wdeclardat_long$idose, wdeclardat_long$DRorDLT2)
wdeclardat_long$label <- NA
wdeclardat_long$label[which(wdeclardat_long$time == max(wdeclardat_long$time ))] <- wdeclardat_long$DRorDLT2[which(wdeclardat_long$time  == max(wdeclardat_long$time ))]


w4 <- wdeclardat_long %>% filter(idose==0.4) %>% 
  ggplot(aes(x=time,y=HazardRate,color=doseandhaz,group=doseandhaz))+
  geom_line(size=2)+
  geom_label_repel(aes(label = label),nudge_x=1,na.rm=TRUE)+
  coord_cartesian(ylim = c(0, 0.055))+
  scale_color_manual(values=c("#9cb99d","#245d67","#003046"))+
  ylab("Hazard Rate")+
  xlab("Weeks")+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),
        panel.background = element_rect(fill = "white"),
        title= element_text("increasing"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  annotate(geom="text", x=40, y=0.05, label="Large Decreasing",
           color="#404040",size=6)


haztraj <- ggarrange(w1, w2, w3, w4)
annotate_figure(haztraj,
                bottom = text_grob("Weeks", color = "#404040", face = "bold", size = 14),
                left = text_grob("Hazard Rate", face="bold",color = "#404040", rot = 90),
                theme(axis.title.y = element_text(margin = margin(t = 0, r = 100, b = 0, l = 0)))
)

#############################################################################







#### Plot Dose-Toxicity
library(tidyverse)
library(readxl)
DTdata <- read_excel("DoseToxicityInfo.xlsx")
DTdata_long <- pivot_longer(DTdata,cols=starts_with("Dose"),names_to = "Dose",values_to = "Prob")
DTdata_long$TypeandTarg <- paste(DTdata_long$ProbType,DTdata_long$TargetDose)

#50/30 Null Const
d1 <- DTdata_long %>% filter(Generator=="Constant" & Scenario=="Null" & (TargetProb==50 | TargetProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Constant: 50/30")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))
  

#50/30 Null Weib Small Dec
d2 <- DTdata_long %>% filter(Generator=="Weibull_sd" & Scenario=="Null" & (TargetProb==50 | TargetProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib SD: 50/30")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#50/30 Null Weib Big Dec
d3 <- DTdata_long %>% filter(Generator=="Weibull_bd" & Scenario=="Null" & (TargetProb==50 | TargetProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib BD: 50/30")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))


#50/30 Null Weib Inc
d4 <- DTdata_long %>% filter(Generator=="Weibull_inc" & Scenario=="Null" & (TargetProb==50 | TargetProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib Inc: 50/30")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#40/20 Null Const
d5 <- DTdata_long %>% filter(Generator=="Constant" & Scenario=="Null" & (TargetProb==40 | TargetProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Constant")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))


#40/20  Null Weib Small Dec
d6 <- DTdata_long %>% filter(Generator=="Weibull_sd" & Scenario=="Null" & (TargetProb==40 | TargetProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib SD")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#40/20  Null Weib Big Dec
d7 <- DTdata_long %>% filter(Generator=="Weibull_bd" & Scenario=="Null" & (TargetProb==40 | TargetProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib LD")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#40/20  Null Weib Inc
d8 <- DTdata_long %>% filter(Generator=="Weibull_inc" & Scenario=="Null" & (TargetProb==40 | TargetProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib Inc")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))


#35/30 Const
d9 <- DTdata_long %>% filter(Generator=="Constant" & Scenario=="LowerDR" & (GenProb==35 | GenProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("35|30 Constant")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))


#35/30 Weib Small Dec
d10 <- DTdata_long %>% filter(Generator=="Weibull_sd" & Scenario=="LowerDR" & (GenProb==35 | GenProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("35|30 Weib SD")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))


#35/30 Weib Big Dec
d11 <- DTdata_long %>% filter(Generator=="Weibull_bd" & Scenario=="LowerDR" & (GenProb==35 | GenProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("35|30 Weib BD")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#35/30 Weib Inc
d12 <- DTdata_long %>% filter(Generator=="Weibull_inc" & Scenario=="LowerDR" & (GenProb==35 | GenProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("35|30 Weib Inc")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#65/30 Const
d13 <- DTdata_long %>% filter(Generator=="Constant" & Scenario=="HigherDR" & (GenProb==65 | GenProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("65|30 Constant")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#65/30 Weib Small Dec
d14 <- DTdata_long %>% filter(Generator=="Weibull_sd" & Scenario=="HigherDR" & (GenProb==65 | GenProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("65|30 Weib SD")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))


#65/30 Weib Big Dec
d15 <- DTdata_long %>% filter(Generator=="Weibull_bd" & Scenario=="HigherDR" & (GenProb==65 | GenProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("65|30 Weib BD")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))


#65/30 Weib inc
d16 <- DTdata_long %>% filter(Generator=="Weibull_inc" & Scenario=="HigherDR" & (GenProb==65 | GenProb==30)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("65|30 Weib Inc")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#25/20 Const
d17 <- DTdata_long %>% filter(Generator=="Constant" & Scenario=="LowerDR" & (GenProb==25 | GenProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Constant: Lower DR")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#25/20 Weib Small Dec
d18 <- DTdata_long %>% filter(Generator=="Weibull_sd" & Scenario=="LowerDR" & (GenProb==25 | GenProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib SD: Lower DR")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))


#25/20 Weib Big Dec
d19 <- DTdata_long %>% filter(Generator=="Weibull_bd" & Scenario=="LowerDR" & (GenProb==25 | GenProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib LD: Lower DR")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))


#25/20 Weib Inc
d20 <- DTdata_long %>% filter(Generator=="Weibull_inc" & Scenario=="LowerDR" & (GenProb==25 | GenProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib Inc: Lower DR")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#55/20 Const
d21 <- DTdata_long %>% filter(Generator=="Constant" & Scenario=="HigherDR" & (GenProb==55 | GenProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Constant: Higher DR")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))


#55/20 Weib Small Dec
d22 <- DTdata_long %>% filter(Generator=="Weibull_sd" & Scenario=="HigherDR" & (GenProb==55 | GenProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib SD: Higher DR")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))



#55/20 Weib Big Dec
d23 <- DTdata_long %>% filter(Generator=="Weibull_bd" & Scenario=="HigherDR" & (GenProb==55 | GenProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Event Probability (%)")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib LD: Higher DR")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))

#55/20 Weib Inc
d24 <- DTdata_long %>% filter(Generator=="Weibull_inc" & Scenario=="HigherDR" & (GenProb==55 | GenProb==20)) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1.5)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=4)+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Dose 5")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14),
        panel.grid.major.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold"))+
  ggtitle("Weib Inc: Higher DR")+
  guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))












DTdata_long$Generator <- factor(DTdata_long$Generator,levels=c("Constant","Weibull_inc","Weibull_sd","Weibull_bd"))

a1 <- DTdata_long %>% filter(GenProb==55 | GenProb==25 | GenProb==20 | GenProb==40) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=1.5)+
  facet_wrap(vars(Generator,Scenario),nrow = 4,scales="free_x")+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Dose 5")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position="none",
        axis.text=element_text(size=8, face="bold"),
        axis.title.y=element_text(hjust = 1,colour="white",face="bold",size=20,margin=margin(0,0,0,5)),
        axis.title.x=element_text(hjust = 1,colour="white",face="bold",size=20,margin=margin(10,0,0,0)),
        panel.grid.major.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold",size=20,margin=margin(0,0,10,0)))+
  #guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))+
  labs(title="",ylab="")
a1


all <- annotate_figure(a1,
                       bottom = text_grob("Dose", vjust=0.6,hjust=-0.5,size=13,face="bold"),
                       right = text_grob(".       .",color="white"),
                       left = text_grob("Toxicity Probability (%)",vjust=1.5,size=13,face="bold", rot = 90))+
  annotate("text", x = 0.25, y = 0.98, label = "Higher DR",size=4,fontface="bold")+
  annotate("text", x = 0.53, y = 0.98, label = "Lower DR",size=4,fontface="bold")+
  annotate("text", x = 0.79, y = 0.98, label = "Matching Dose",size=4,fontface="bold")+
  annotate("text", x=0.96,y=0.85,label="Constant",size=4,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.62,label="Small\nIncreasing",size=4,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.41,label="Small\nDecreasing",size=4,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.18,label="Large\nDecreasing",size=4,color="#292929",angle="270")+
  annotate("text", x=0.84,y=0.01,label="DLT              DR",size=3,color="#292929")+
  annotate("text", x=0.80,y=0.02,label="__      _ _",size=6,fontface="bold",color="#003046")+
  annotate("text", x=0.145,y=0.03,label="Target DLT",size=3,color="#292929")+
  annotate("text", x=0.24,y=0.01,label="Dose 1        Dose 3        Dose 5",size=3,color="#292929")+
  annotate("text", x=0.115,y=0.033,label=".",size=20,fontface="bold",color="#003046")+
  annotate("text", x=0.20,y=0.033,label=".",size=20,fontface="bold",color="#9cb99d")+
  annotate("text", x=0.285,y=0.033,label=".",size=20,fontface="bold",color="#B5A573")
all
pdf("dosetoxtraj4020.pdf", width=8, height=10); all; dev.off()



















a2 <- DTdata_long %>% filter(GenProb==65 | GenProb==35 | GenProb==30 | GenProb==50) %>% 
  ggplot(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg))+
  geom_line(aes(linetype=ProbType),size=1)+
  geom_point(aes(x=Dose,y=Prob,color=TypeandTarg,group=TypeandTarg),size=1.5)+
  facet_wrap(vars(Generator,Scenario),nrow = 4,scales="free_x")+
  scale_color_manual(values=c("#9cb99d","#B5A573","#003046","#9cb99d","#B5A573","#003046"))+
  coord_cartesian(ylim = c(0, 100))+
  ylab("Dose 5")+
  xlab(element_blank())+
  theme(panel.grid.major.y = element_line(colour="#404040",linetype=2),panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position="none",
        axis.text=element_text(size=8, face="bold"),
        axis.title.y=element_text(hjust = 1,colour="white",face="bold",size=20,margin=margin(0,0,0,5)),
        axis.title.x=element_text(hjust = 1,colour="white",face="bold",size=20,margin=margin(10,0,0,0)),
        panel.grid.major.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 1,colour="#404040",face="bold",size=20,margin=margin(0,0,10,0)))+
  #guides(col = "none")+
  scale_x_discrete(labels=c("1","2","3","4","5"))+
  labs(title="",ylab="")
a2


all2 <- annotate_figure(a2,
                       bottom = text_grob("Dose", vjust=0.6,hjust=-0.5,size=13,face="bold"),
                       right = text_grob(".       .",color="white"),
                       left = text_grob("Toxicity Probability (%)",vjust=1.5,size=13,face="bold", rot = 90))+
  annotate("text", x = 0.25, y = 0.98, label = "Higher DR",size=4,fontface="bold")+
  annotate("text", x = 0.53, y = 0.98, label = "Lower DR",size=4,fontface="bold")+
  annotate("text", x = 0.79, y = 0.98, label = "Matching Dose",size=4,fontface="bold")+
  annotate("text", x=0.96,y=0.85,label="Constant",size=4,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.62,label="Small\nIncreasing",size=4,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.41,label="Small\nDecreasing",size=4,color="#292929",angle="270")+
  annotate("text", x=0.96,y=0.18,label="Large\nDecreasing",size=4,color="#292929",angle="270")+
  annotate("text", x=0.84,y=0.01,label="DLT              DR",size=3,color="#292929")+
  annotate("text", x=0.80,y=0.02,label="__      _ _",size=6,fontface="bold",color="#003046")+
  annotate("text", x=0.145,y=0.03,label="Target DLT",size=3,color="#292929")+
  annotate("text", x=0.24,y=0.01,label="Dose 1        Dose 3        Dose 5",size=3,color="#292929")+
  annotate("text", x=0.115,y=0.033,label=".",size=20,fontface="bold",color="#003046")+
  annotate("text", x=0.20,y=0.033,label=".",size=20,fontface="bold",color="#9cb99d")+
  annotate("text", x=0.285,y=0.033,label=".",size=20,fontface="bold",color="#B5A573")
all2
pdf("dosetoxtraj.pdf", width=8, height=10); all2; dev.off()
