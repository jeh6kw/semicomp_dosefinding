library(tidyverse)
library(readxl)
library(gridExtra)
library(ggpubr)
setwd("C:/Users/jesse/Documents/Coursework/Research/Conaway_Fall22")
####################################
############ Efficacy Estimand Examples
####################################
set.seed(100)



### Big 
bigdat <- read_excel("efficacydataexamples.xlsx"); bigdat$Patient <- as.character(bigdat$Patient)
bigdatcens <- bigdat %>% filter(CENS=="Yes")
bigdatDLT <-  bigdat %>% filter(DLT=="Yes")
bigdatDR <-  bigdat %>% filter(DR=="Yes")
bigdat$Dose <- factor(bigdat$Dose,levels=c("one","two","three"))

bg1 <- ggplot(bigdat, aes(x = Time, y = CfB,group = Patient, color=Dose)) + 
  geom_line(linewidth = .1)+
  geom_point(size=0.5)+
  geom_hline(yintercept=0,linetype=1,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-30,linetype=2,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-80,linetype=2,color="darkgrey",linewidth=0.1)+
  coord_cartesian(ylim=c(-100, 70),xlim=c(2.3,52))+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        axis.line = element_line(colour ="darkgrey",linewidth=0.1),
        legend.position = "right",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5,color='#3A3838',size=10),
        legend.title = element_text(color="#3A3838",size=8),
        legend.text = element_text(color="#3A3838",size=8),
        legend.key=element_blank())+
  scale_color_manual(values=c('#6C986D','#ff6361','#002332'),name="Initial Dose")+
  ggtitle("Large Regression")+
  geom_point(data=bigdatcens,shape=19,size=2)+
  geom_point(data=bigdatDLT,shape=15,size=2)+
  geom_point(data=bigdatDR,shape=17,size=1.5)
bg1



### Medium 
meddat <-  read_excel("efficacydataexamples.xlsx",sheet=2)
meddatcens <- meddat %>% filter(CENS=="Yes")
meddatDLT <-  meddat %>% filter(DLT=="Yes")
meddatDR <-  meddat %>% filter(DR=="Yes")
meddat$Dose <- factor(meddat$Dose,levels=c("one","two","three"))

mg1 <- ggplot(meddat, aes(x = Time, y = CfB,group = Patient, color=Dose)) + 
  geom_line(linewidth=0.1)+
  geom_point(size=0.5)+
  geom_hline(yintercept=0,linetype=1,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-30,linetype=2,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-80,linetype=2,color="darkgrey",linewidth=0.1)+
  coord_cartesian(ylim=c(-100, 70),xlim=c(2.3,52))+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        axis.line = element_line(colour ="darkgrey",linewidth=0.1),
        legend.position = "right",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(color="#3A3838",size=10,margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5,color='#3A3838',size=10),
        legend.title = element_text(color="#3A3838",size=8),
        legend.text = element_text(color="#3A3838",size=8),
        legend.key=element_blank())+
  scale_color_manual(values=c('#6C986D','#ff6361','#002332'),name="Initial Dose")+
  ggtitle("Medium Regression")+
  geom_point(data=meddatcens,shape=19,size=2)+
  geom_point(data=meddatDLT,shape=15,size=2)+
  geom_point(data=meddatDR,shape=17,size=1.5)+
  labs(y = "")
mg1





### Small
smalldat <-  read_excel("efficacydataexamples.xlsx",sheet=3)
smalldatcens <- smalldat %>% filter(CENS=="Yes")
smalldatDLT <-  smalldat %>% filter(DLT=="Yes")
smalldatDR <-  smalldat %>% filter(DR=="Yes")
smalldat$Dose <- factor(smalldat$Dose,levels=c("one","two","three"))


sg1 <- ggplot(smalldat, aes(x = Time, y = CfB,group = Patient, color=Dose)) + 
  geom_line(linewidth=0.1)+
  geom_point(size=.5)+
  geom_hline(yintercept=0,linetype=1,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-30,linetype=2,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-80,linetype=2,color="darkgrey",linewidth=0.1)+
  coord_cartesian(ylim=c(-100, 70),xlim=c(2.3,52))+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        axis.line = element_line(colour ="darkgrey",linewidth=0.1),
        legend.position = "right",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5,color='#3A3838',size=10),
        legend.title = element_text(color="#3A3838",size=8),
        legend.text = element_text(color="#3A3838",size=8),
        legend.key=element_blank())+
  scale_color_manual(values=c('#6C986D','#ff6361','#002332'),name="Initial Dose")+
  ggtitle("Small Regression")+
  geom_point(data=smalldatcens,shape=19,size=2)+
  geom_point(data=smalldatDLT,shape=15,size=2)+
  geom_point(data=smalldatDR,shape=17,size=1.5)
sg1








##############  IMPUTE DATA  ###############

#Small Noise
library(mice)
library(miceadds)
library(micemd)
smallmissdat <-  read_excel("efficacydataexamples.xlsx",sheet=6)
smallmissdat$Time2 <- smallmissdat$Time^2
smallmissdat$ImpDose <- as.factor(smallmissdat$ImpDose)
smallmissdat$Dose <- as.factor(smallmissdat$Dose)


imp <- mice(smallmissdat)
impmat <- imp$predictorMatrix 
impmat[,c("Dose")] <- 0;
impmat[,c("Patient")] <- -2; impmat[,c("Time")] <- 0;impmat[,c("Time2")] <- 0;
impmat[,c("CfB")] <- c(1,1,2,0,1,2); impmat[,c("ImpDose")] <- c(1,1,2,1,0,2)
impmat

methimp <- imp$method
methimp["CfB"] <- "2l.pan"
methimp

micelongs <- mice(smallmissdat, meth = methimp, pred = impmat, m= 20,
                 maxit = 20, seed = 123, printFlag = FALSE)

simpscfb <- apply(micelongs$imp$CfB,1,function(x) median(x))
simpsimpdose <- apply(micelongs$imp$ImpDose,1,function(x) names(which.max(table(x))))


#Medium Noise
set.seed(123)
library(mice)
medmissdat <-  read_excel("efficacydataexamples.xlsx",sheet=5)
medmissdat$Time2 <- medmissdat$Time^2
medmissdat$ImpDose <- as.factor(medmissdat$ImpDose)
medmissdat$Dose <- as.factor(medmissdat$Dose)


imp <- mice(medmissdat)
impmat <- imp$predictorMatrix 
impmat[,c("Dose")] <- 0; impmat[,c("Patient")] <- -2; impmat[,c("Time")] <- 0; impmat[,c("Time2")] <- 0;
impmat[,c("CfB")] <- c(1,1,2,0,1,2); impmat[,c("ImpDose")] <- c(1,1,2,1,0,2)
impmat

methimp <- imp$method
methimp["CfB"] <- "2l.pan"
methimp

micelongm <- mice(medmissdat, meth = methimp, pred = impmat, m= 20,
                 maxit = 20, seed = 123, printFlag = FALSE)


mimpscfb <- apply(micelongm$imp$CfB,1,function(x) median(x))
mimpsimpdose <- apply(micelongm$imp$ImpDose,1,function(x) names(which.max(table(x))))



#Big Noise
library(mice)
bigmissdat <- read_excel("efficacydataexamples.xlsx",sheet=4)
bigmissdat$Time2 <- bigmissdat$Time^2
bigmissdat$ImpDose <- as.factor(bigmissdat$ImpDose)
bigmissdat$Dose <- as.factor(bigmissdat$Dose)


imp <- mice(bigmissdat)
impmat <- imp$predictorMatrix 
impmat[,c("Dose")] <- 0
impmat[,c("Patient")] <- -2
impmat[,c("Time")] <- 0
impmat[,c("Time2")] <- 0
impmat[,c("CfB")] <- c(1,1,2,0,1,2); impmat[,c("ImpDose")] <- c(1,1,2,1,0,2)
impmat

methimp <- imp$method
methimp["CfB"] <- "2l.pan"
methimp

micelongb <- mice(bigmissdat, meth = methimp, pred = impmat, m= 20,
                 maxit = 20, seed = 123, printFlag = FALSE)


bimpscfb <- apply(micelongb$imp$CfB,1,function(x) median(x))
bimpsimpdose <- apply(micelongb$imp$ImpDose,1,function(x) names(which.max(table(x))))








############################################################################
############################################################################



### Small Imp
smalldatimp <- smallmissdat
smalldatimp$CfB[as.numeric(names(simpscfb))] <- simpscfb
smalldatimp$ImpDose[as.numeric(names(simpsimpdose))] <- simpsimpdose
#smalldatDR <-  smalldatimp %>% filter(DRorDLT=="Yes")
smalldatimp$Dose <- factor(smalldatimp$Dose,levels=c("one","two","three"))

sg1
sg1_imp <- ggplot(smalldatimp, aes(x = Time, y = CfB,group = Patient, color=Dose)) + 
  geom_line(linetype = "dashed",linewidth=0.1)+
  geom_line(data=smalldat,aes(x=Time,y=CfB,group=Patient,color=Dose),size=0.2)+
  geom_point(size=.5)+
  geom_hline(yintercept=0,linetype=1,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-30,linetype=2,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-80,linetype=2,color="darkgrey",linewidth=0.1)+
  coord_cartesian(ylim=c(-100, 70),xlim=c(2.3,52))+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        axis.line = element_line(colour ="darkgrey",linewidth=0.1),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5,color="#3A3838",size=10),
        legend.title = element_text(color="#3A3838",size=8),
        legend.text = element_text(color="#3A3838",size=8),
        legend.key=element_blank())+
  scale_color_manual(values=c('#6C986D','#ff6361','#002332'),name="Initial Dose")+
  ggtitle("Small Regression")+
  geom_point(data=smalldatcens,shape=19,size=2)+
  geom_point(data=smalldatDLT,shape=15,size=2)+
  geom_point(data=smalldatDR,shape=17,size=1.5)
sg1_imp


### Medium Imp
meddatimp <- medmissdat
meddatimp$CfB[as.numeric(names(mimpscfb))] <- mimpscfb
meddatimp$ImpDose[as.numeric(names(mimpsimpdose))] <- mimpsimpdose
#meddatDR <-  meddatimp %>% filter(DRorDLT=="Yes")
meddatimp$Dose <- factor(meddatimp$Dose,levels=c("one","two","three"))

mg1
mg1_imp <- ggplot(meddatimp, aes(x = Time, y = CfB,group = Patient, color=Dose)) + 
  geom_line(linetype = "dashed",linewidth=0.1)+
  geom_line(data=meddat,aes(x=Time,y=CfB,group=Patient,color=Dose),size=0.2)+
  geom_point(size=.5)+
  geom_hline(yintercept=0,linetype=1,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-30,linetype=2,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-80,linetype=2,color="darkgrey",linewidth=0.1)+
  coord_cartesian(ylim=c(-100, 70),xlim=c(2.3,52))+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        axis.line = element_line(colour ="darkgrey",linewidth=0.1),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(color="#3A3838",margin = margin(t = 0, r = 15, b = 0, l = 0),size=10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5,color='#3A3838',size=10),
        legend.title = element_text(color="#3A3838",size=8),
        legend.text = element_text(color="#3A3838",size=8),
        legend.key=element_blank())+
  scale_color_manual(values=c('#6C986D', '#ff6361','#002332'),name="Initial Dose")+
  ggtitle("Medium Regression")+
  geom_point(data=meddatcens,shape=19,size=2)+
  geom_point(data=meddatDLT,shape=15,size=2)+
  geom_point(data=meddatDR,shape=17,size=1.5)+
  labs(y="")
mg1_imp


### Big Imp
bigdatimp <- bigmissdat
bigdatimp$CfB[as.numeric(names(bimpscfb))] <- bimpscfb
bigdatimp$ImpDose[as.numeric(names(bimpsimpdose))] <- bimpsimpdose
#bigdatDR <-  bigdatimp %>% filter(DRorDLT=="Yes")
bigdatimp$Dose <- factor(bigdatimp$Dose,levels=c("one","two","three"))

bg1
bg1_imp <- ggplot(bigdatimp, aes(x = Time, y = CfB,group = Patient, color=Dose)) + 
  geom_line(linetype = "dashed",linewidth=0.1)+
  geom_line(data=bigdat,aes(x=Time,y=CfB,group=Patient,color=Dose),size=.2)+
  geom_point(size=0.5)+
  geom_hline(yintercept=0,linetype=1,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-30,linetype=2,color="darkgrey",linewidth=0.1)+
  geom_hline(yintercept=-80,linetype=2,color="darkgrey",linewidth=0.1)+
  coord_cartesian(ylim=c(-100, 70),xlim=c(2.3,52))+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#F4F5F0"),
        axis.line = element_line(colour ="darkgrey",linewidth=0.1),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5,color="#3A3838",size=10),
        legend.title = element_text(color="#3A3838",size=8),
        legend.text = element_text(color="#3A3838",size=8),
        legend.key=element_blank())+
  scale_color_manual(values=c('#6C986D', '#ff6361','#002332'),name="Initial Dose")+
  ggtitle("Large Regression")+
  geom_point(data=bigdatcens,shape=19,size=2)+
  geom_point(data=bigdatDLT,shape=15,size=2)+
  geom_point(data=bigdatDR,shape=17,size=1.5)
bg1_imp



































#############################################################################################
##############################################################################################
#######################            PLOT ORGANIZATION   #######################################

library(ggpubr)


p1 <- ggarrange(sg1 ,mg1, bg1,
                labels = NULL,
                ncol = 1,
                align = "v", 
                common.legend = TRUE,
                legend="right")
p1 <- annotate_figure(p1, 
                      bottom = text_grob("Time (weeks)", vjust=0.5,hjust=0.55,size=10,color="#3A3838",face="bold"),
                      left = text_grob("Change in tumor burden from baseline (%)", vjust=2,hjust=0.55,size=10,color="#3A3838",face="bold",rot=90))
p1


pdf("nonimp_efftraj.pdf", width=4, height=8); p1; dev.off()






p2 <- ggarrange(sg1_imp ,mg1_imp, bg1_imp,
                labels = NULL,
                ncol = 1,
                align = "v", 
                common.legend = TRUE,
                legend="right")
p2 <- annotate_figure(p2, 
                      bottom = text_grob("Time (weeks)", vjust=0.5,hjust=0.55,size=10,color="#3A3838",face="bold"),
                      left = text_grob("Change in tumor burden from baseline (%)", vjust=2,hjust=0.55,size=10,face="bold",color="#3A3838",rot=90))
p2


pdf("imp_efftraj.pdf", width=4, height=8)
p2
dev.off()
