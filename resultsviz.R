setwd("C:/Users/jesse/Documents/Coursework/Research/Conaway_Fall22/Rivannaresults_tox")
library(readxl)
res <- read_excel("SummaryResults_AllMethods.xlsx")
library(tidyverse)
library(ggplot2)
res_red <- res %>% 
  filter(dlt_target == 20) %>% 
  #filter(dlt_target == 30) %>% 
  #filter((trueDLT_dose == 1 & DRsituation == "null") |
           # (trueDLT_dose == 3) |
           # (trueDLT_dose == 5 & DRsituation == "null") |
           # (trueDLT_dose == 5 & DRsituation =="higherDR")) %>%
  filter(method == "constant_uninform")  
  







############################# PCS and Overdose

res_po <- res_red %>% 
  select(DRsituation, datagen, trueDLT_dose, n, PCS_myconst, PCS_myconst_cutoff,	PCS_myconst_31,
         PCS_myweib,	PCS_myweib_cutoff,	PCS_myweib_cutoff31,
         overdose_myconst,	overdose_myconst_cutoff,	overdose_myconst_31,	
         overdose_myweib,	overdose_myweib_cutoff, overdose_myweib_cutoff31)
res_polong <- res_po %>% 
  pivot_longer(PCS_myconst:overdose_myweib_cutoff31,names_to="model",values_to="value") %>% 
  separate(model,into=c("measure","modeltype"),sep="_",extra="merge")

comp <- res_polong %>% group_by(DRsituation,datagen, trueDLT_dose,n ,modeltype) %>% summarize(value = sum(value))
comp$measure <- "fill"; comp$value <- 100-comp$value
res_polong <- rbind(res_polong,comp)
res_polong$measure <- factor(res_polong$measure,levels=c("overdose","fill","PCS"))
res_polong$datagen <- factor(res_polong$datagen,levels=c("constant","weib_inc","weib_dec_small","weib_dec_big"))

##n = 30, scenario 1
n30_1 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 1")

##n = 30, scenario 2
n30_2 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 2")


##n = 30, scenario 3
n30_3 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 3")

##n = 30, scenario 5
n30_5 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 5")

##n = 30, scenario 8
n30_8 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 8")

##n = 30, scenario 9
n30_9 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 9")+
  labs(y=NULL)


library(ggpubr)
library(grid)
s2 <- ggarrange(n30_1,n30_2,n30_3,n30_5,n30_8,n30_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
s2

my2 <- annotate_figure(s2,
                        bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.22,y=0.01,label="PCS        Overdose",size=3,color="#292929")+
  annotate("text", x=0.13,y=0.033,label=".",size=20,fontface="bold",color="#9cb99d")+
  annotate("text", x=0.207,y=0.033,label=".",size=20,fontface="bold",color="#003046")
my2

pdf("ex1.pdf", width=2, height=3); my2; dev.off()


















##n = 60, scenario 1
n60_1 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 1")

##n = 30, scenario 2
n60_2 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 2")


##n = 60, scenario 3
n60_3 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 3")

##n = 60, scenario 5
n60_5 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 5")

##n = 60, scenario 8
n60_8 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 8")

##n = 60, scenario 9
n60_9 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
s3 <- ggarrange(n60_1,n60_2,n60_3,n60_5,n60_8,n60_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
s3

my3 <- annotate_figure(s3,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.22,y=0.01,label="PCS        Overdose",size=3,color="#292929")+
  annotate("text", x=0.13,y=0.033,label=".",size=20,fontface="bold",color="#9cb99d")+
  annotate("text", x=0.207,y=0.033,label=".",size=20,fontface="bold",color="#003046")
my3

pdf("pcsover_sec1_n60.pdf", width=7, height=10); my3; dev.off()




















##n = 100, scenario 1
n100_1 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 1")

##n = 100, scenario 2
n100_2 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 2")


##n = 100, scenario 3
n100_3 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 3")

##n = 100, scenario 5
n100_5 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 5")

##n = 100, scenario 8
n100_8 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 8")

##n = 100, scenario 9
n100_9 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",width = 0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
s4 <- ggarrange(n100_1,n100_2,n100_3,n100_5,n100_8,n100_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
s4

my4 <- annotate_figure(s4,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.22,y=0.01,label="PCS        Overdose",size=3,color="#292929")+
  annotate("text", x=0.13,y=0.033,label=".",size=20,fontface="bold",color="#9cb99d")+
  annotate("text", x=0.207,y=0.033,label=".",size=20,fontface="bold",color="#003046")
my4

pdf("pcsover_sec1_n100.pdf", width=7, height=10); my4; dev.off()















#DRs are gold
################################### Acc Indices
res_po <- res_red %>% 
  select(DRsituation, datagen,trueDLT_dose, n, AccIndexDLT_myconst,	AccIndexDLT_myconst_cutoff,	
         AccIndexDLT_myconst_31,	AccIndexDLT_myweib, AccIndexDLT_myweib_cutoff,AccIndexDLT_myweib_cutoff31,
         AccIndexDR_myconst,	AccIndexDR_myconst_cutoff,	AccIndexDR_myconst_31,
         AccIndexDR_myweib,	AccIndexDR_myweib_cutoff,	AccIndexDR_myweib_cutoff31)
res_polong <- res_po %>% 
  pivot_longer(AccIndexDLT_myconst:AccIndexDR_myweib_cutoff31,names_to="model",values_to="value") %>% 
  separate(model,into=c("measure","modeltype"),sep="_",extra="merge")


res_polong$measure <- factor(res_polong$measure,levels=c("AccIndexDR","fill","AccIndexDLT"))
res_polong$datagen <- factor(res_polong$datagen,levels=c("constant","weib_inc","weib_dec_small","weib_dec_big"))



##n = 30, scenario 1
z30_1 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 1")


##n = 30, scenario 2
z30_2 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==3 )%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 2")


##n = 30, scenario 3
z30_3 <-  res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 3")

##n = 30, scenario 5
z30_5 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3 & measure=="AccIndexDLT")%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 5")

##n = 30, scenario 8
z30_8 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 8")

##n = 30, scenario 9
z30_9 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
z2 <- ggarrange(z30_1,z30_2,z30_3,z30_5,z30_8,z30_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
z2

fe2 <- annotate_figure(z2,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.33,y=0.01,label="Accuracy Index: DR        Accuracy Index: DLT",size=3,color="#292929")+
  annotate("text", x=0.13,y=0.033,label=".",size=20,fontface="bold",color="#c6c88c")+
  annotate("text", x=0.33,y=0.033,label=".",size=20,fontface="bold",color="#736A35")
fe2

pdf("drdlt_sec1_n30.pdf", width=7, height=10); fe2; dev.off()












##n = 60, scenario 1
z60_1 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 1")


##n = 60, scenario 2
z60_2 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==3 )%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 2")


##n = 60, scenario 3
z60_3 <-  res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 3")

##n = 60, scenario 5
z60_5 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3 & measure=="AccIndexDLT")%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 5")

##n = 60, scenario 8
z60_8 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 8")

##n = 60, scenario 9
z60_9 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
z3 <- ggarrange(z60_1,z60_2,z60_3,z60_5,z60_8,z60_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
z3

fe3 <- annotate_figure(z3,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.33,y=0.01,label="Accuracy Index: DR        Accuracy Index: DLT",size=3,color="#292929")+
  annotate("text", x=0.13,y=0.033,label=".",size=20,fontface="bold",color="#c6c88c")+
  annotate("text", x=0.33,y=0.033,label=".",size=20,fontface="bold",color="#736A35")
fe3

pdf("drdlt_sec1_n60.pdf", width=7, height=10); fe3; dev.off()








##n = 100, scenario 1
z100_1 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 1")


##n = 100, scenario 2
z100_2 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==3 )%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 2")


##n = 100, scenario 3
z100_3 <-  res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 3")

##n = 100, scenario 5
z100_5 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3 & measure=="AccIndexDLT")%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 5")

##n = 100, scenario 8
z100_8 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 8")

##n = 100, scenario 9
z100_9 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=measure))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="#585858"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("const","const 3:1","const/mutd","weib","weib/mutd","weib 3:1/mutd"))+
  scale_fill_manual(values=c('#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
z4 <- ggarrange(z100_1,z100_2,z100_3,z100_5,z100_8,z100_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
z4

fe4 <- annotate_figure(z4,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.33,y=0.01,label="Accuracy Index: DR        Accuracy Index: DLT",size=3,color="#292929")+
  annotate("text", x=0.13,y=0.033,label=".",size=20,fontface="bold",color="#c6c88c")+
  annotate("text", x=0.33,y=0.033,label=".",size=20,fontface="bold",color="#736A35")
fe4

pdf("drdlt_sec1_n100.pdf", width=7, height=10); fe4; dev.off()








































































#####################################################################################################







############################# PCS and Overdose, Comparator Models

res_po <- res_red %>% 
  select(DRsituation, datagen, trueDLT_dose, n, PCS_myweib_cutoff, PCS_titecrm, PCS_titecrmmcext, PCS_crmmcext,
         overdose_myweib_cutoff, overdose_titecrm, overdose_titecrmmcext,overdose_crmmcext)
res_polong <- res_po %>% 
  pivot_longer(PCS_myweib_cutoff:overdose_crmmcext,names_to="model",values_to="value") %>% 
  separate(model,into=c("measure","modeltype"),sep="_",extra="merge")

comp <- res_polong %>% group_by(DRsituation,datagen, trueDLT_dose,n ,modeltype) %>% summarize(value = sum(value))
comp$measure <- "fill"; comp$value <- 100-comp$value
res_polong <- rbind(res_polong,comp)
res_polong$measure <- factor(res_polong$measure,levels=c("overdose","fill","PCS"))
res_polong$datagen <- factor(res_polong$datagen,levels=c("constant","weib_inc","weib_dec_small","weib_dec_big"))
res_polong$modeltype <- factor(res_polong$modeltype,levels=c("titecrm","crmmcext","titecrmmcext","myweib_cutoff"))
res_polong$colorcombo <- paste0(res_polong$measure,res_polong$modeltype)
res_polong$colorcombo <- factor(res_polong$colorcombo,levels = c("overdosetitecrm","filltitecrm","PCStitecrm","overdosecrmmcext","fillcrmmcext","PCScrmmcext","overdosetitecrmmcext","filltitecrmmcext","PCStitecrmmcext","overdosemyweib_cutoff","fillmyweib_cutoff","PCSmyweib_cutoff"))

##n = 30, scenario 1
c30_1 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 1")



##n = 30, scenario 2
c30_2 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 2")


##n = 30, scenario 3
c30_3 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 3")

##n = 30, scenario 5
c30_5 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 5")

##n = 30, scenario 8
c30_8 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 8")

##n = 30, scenario 9
c30_9 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
c2 <- ggarrange(c30_1,c30_2,c30_3,c30_5,c30_8,c30_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
c2

it2 <- annotate_figure(c2,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.345,y=0.01,label="PCS        PCS (Comparator Model)        Overdose",size=3,color="#292929")+
  annotate("text", x=0.135,y=0.033,label=".",size=20,fontface="bold",color="#9cb99d")+
  annotate("text", x=0.21,y=0.033,label=".",size=20,fontface="bold",color="#D0DED0")+
  annotate("text", x=0.45,y=0.033,label=".",size=20,fontface="bold",color="#003046")
it2

pdf("pcsover_2_n30.pdf", width=7, height=10); it2; dev.off()
















##n = 60, scenario 1
c60_1 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 1")



##n = 60, scenario 2
c60_2 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 2")


##n = 60, scenario 3
c60_3 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 3")

##n = 60, scenario 5
c60_5 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 5")

##n = 60, scenario 8
c60_8 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 8")

##n = 60, scenario 9
c60_9 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
c3 <- ggarrange(c60_1,c60_2,c60_3,c60_5,c60_8,c60_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
c3

it3 <- annotate_figure(c3,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.345,y=0.01,label="PCS        PCS (Comparator Model)        Overdose",size=3,color="#292929")+
  annotate("text", x=0.135,y=0.033,label=".",size=20,fontface="bold",color="#9cb99d")+
  annotate("text", x=0.21,y=0.033,label=".",size=20,fontface="bold",color="#D0DED0")+
  annotate("text", x=0.45,y=0.033,label=".",size=20,fontface="bold",color="#003046")
it3

pdf("pcsover_2_n60.pdf", width=7, height=10); it3; dev.off()



















##n = 100, scenario 1
c100_1 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 1")



##n = 100, scenario 2
c100_2 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 2")


##n = 100, scenario 3
c100_3 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 3")

##n = 100, scenario 5
c100_5 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 5")

##n = 100, scenario 8
c100_8 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 8")

##n = 100, scenario 9
c100_9 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046','#F2F2F2','#D0DED0','#003046', '#F2F2F2', '#9cb99d'))+
  geom_hline(yintercept=c(25,50,75),linetype=3,color="#585858")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
c4 <- ggarrange(c100_1,c100_2,c100_3,c100_5,c100_8,c100_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
c4

it4 <- annotate_figure(c4,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.345,y=0.01,label="PCS        PCS (Comparator Model)        Overdose",size=3,color="#292929")+
  annotate("text", x=0.135,y=0.033,label=".",size=20,fontface="bold",color="#9cb99d")+
  annotate("text", x=0.21,y=0.033,label=".",size=20,fontface="bold",color="#D0DED0")+
  annotate("text", x=0.45,y=0.033,label=".",size=20,fontface="bold",color="#003046")
it4

pdf("pcsover_2_n100.pdf", width=7, height=10); it4; dev.off()

















#DRs are gold
################################### Acc Indices, model comparators
res_po <- res_red %>% 
  select(DRsituation, datagen,trueDLT_dose, n, AccIndexDLT_myweib_cutoff, AccIndexDLT_titecrm, AccIndexDLT_titecrmmcext,
         AccIndexDLT_crmmcext, AccIndexDR_myweib_cutoff, AccIndexDR_titecrm, AccIndexDR_titecrmmcext, AccIndexDR_crmmcext)
res_polong <- res_po %>% 
  pivot_longer(AccIndexDLT_myweib_cutoff:AccIndexDR_crmmcext,names_to="model",values_to="value") %>% 
  separate(model,into=c("measure","modeltype"),sep="_",extra="merge")


res_polong$measure <- factor(res_polong$measure,levels=c("AccIndexDR","fill","AccIndexDLT"))
res_polong$datagen <- factor(res_polong$datagen,levels=c("constant","weib_inc","weib_dec_small","weib_dec_big"))
res_polong$modeltype <- factor(res_polong$modeltype, levels = c("titecrm","crmmcext","titecrmmcext","myweib_cutoff"))
res_polong$colorcombo <- paste0(res_polong$measure,res_polong$modeltype)
res_polong$colorcombo <- factor(res_polong$colorcombo,levels = c("AccIndexDRtitecrm","AccIndexDLTtitecrm","AccIndexDRcrmmcext","AccIndexDLTcrmmcext",
                                                                 "AccIndexDRtitecrmmcext","AccIndexDLTtitecrmmcext","AccIndexDRmyweib_cutoff","AccIndexDLTmyweib_cutoff"))


##n = 30, scenario 1
o30_1 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE','#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 1")


##n = 30, scenario 2
o30_2 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==3 )%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE','#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 2")


##n = 30, scenario 3
o30_3 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE','#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 3")

##n = 30, scenario 5
o30_5 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3 & measure=="AccIndexDLT")%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#D0CECE', '#D0CECE', '#D0CECE','#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 5")

##n = 30, scenario 8
o30_8 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6","#E7E6E6","#E7E6E6",'#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 8")

##n = 30, scenario 9
o30_9 <- res_polong %>% filter(n==30) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #cale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6","#E7E6E6","#E7E6E6",'#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
o2 <- ggarrange(o30_1,o30_2,o30_3,o30_5,o30_8,o30_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
o2

mo2 <- annotate_figure(o2,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.36,y=0.01,label="Accuracy Index: DR              Accuracy Index: DLT",size=3,color="#292929")+
  annotate("text", x=0.13,y=0.033,label=".",size=20,fontface="bold",color="#c6c88c")+
  annotate("text", x=0.15,y=0.033,label=".",size=20,fontface="bold",color="#E7E6E6")+
  annotate("text", x=0.36,y=0.033,label=".",size=20,fontface="bold",color="#736A35")+
  annotate("text", x=0.38,y=0.033,label=".",size=20,fontface="bold",color="#D0CECE")
mo2

pdf("drdlt_sec2_n30.pdf", width=7, height=10); mo2; dev.off()










##n = 60, scenario 1
o60_1 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE','#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 1")


##n = 60, scenario 2
o60_2 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==3 )%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE','#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 2")


##n = 60, scenario 3
o60_3 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE','#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 3")

##n = 60, scenario 5
o60_5 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3 & measure=="AccIndexDLT")%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#D0CECE', '#D0CECE', '#D0CECE','#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 5")

##n = 60, scenario 8
o60_8 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6","#E7E6E6","#E7E6E6",'#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 8")

##n = 60, scenario 9
o60_9 <- res_polong %>% filter(n==60) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #cale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6","#E7E6E6","#E7E6E6",'#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
o3 <- ggarrange(o60_1,o60_2,o60_3,o60_5,o60_8,o60_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
o3

mo3 <- annotate_figure(o3,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.36,y=0.01,label="Accuracy Index: DR              Accuracy Index: DLT",size=3,color="#292929")+
  annotate("text", x=0.13,y=0.033,label=".",size=20,fontface="bold",color="#c6c88c")+
  annotate("text", x=0.15,y=0.033,label=".",size=20,fontface="bold",color="#E7E6E6")+
  annotate("text", x=0.36,y=0.033,label=".",size=20,fontface="bold",color="#736A35")+
  annotate("text", x=0.38,y=0.033,label=".",size=20,fontface="bold",color="#D0CECE")
mo3

pdf("drdlt_sec2_n60.pdf", width=7, height=10); mo3; dev.off()








##n = 100, scenario 1
o100_1 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==1) %>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE','#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 1")


##n = 100, scenario 2
o100_2 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==3 )%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE','#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 2")


##n = 100, scenario 3
o100_3 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="null" & trueDLT_dose==5)%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE',"#E7E6E6",'#D0CECE','#c6c88c', '#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 3")

##n = 100, scenario 5
o100_5 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="lowerDR" & trueDLT_dose==3 & measure=="AccIndexDLT")%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c('#D0CECE', '#D0CECE', '#D0CECE','#736A35'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 5")

##n = 100, scenario 8
o100_8 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="higherDR" & trueDLT_dose==3 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right")+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        strip.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  scale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6","#E7E6E6","#E7E6E6",'#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 8")

##n = 100 scenario 9
o100_9 <- res_polong %>% filter(n==100) %>% filter(DRsituation=="higherDR" & trueDLT_dose==5 & measure=="AccIndexDR")%>% 
  ggplot(aes(x=modeltype,y=value,fill=colorcombo))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~datagen,
             nrow=4,
             strip.position = "right",
             labeller = as_labeller(c("constant" = "Constant","weib_dec_big" = "Large Decreasing","weib_dec_small" = "Small Decreasing","weib_inc" = "Small Increasing"),
                                    default = label_wrap_gen(width=1)))+
  coord_flip()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        strip.background  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(colour = '#585858',size=8),
        plot.title = element_text(hjust = 0.5,face="bold",margin=margin(5,0,5,0)))+
  #cale_x_discrete(labels=c("TITECRM","CRM-MC","TITECRM-MC","W+M"))+
  scale_fill_manual(values=c("#E7E6E6","#E7E6E6","#E7E6E6",'#c6c88c'))+
  geom_hline(yintercept=c(-1,-0.75,-0.5,-0.25,0.25,.50,.75),linetype=3,color="grey")+
  ggtitle("Scenario 9")


library(ggpubr)
library(grid)
o4 <- ggarrange(o100_1,o100_2,o100_3,o100_5,o100_8,o100_9,
                labels = NULL,
                ncol = 2,
                nrow=3)
o4

mo4 <- annotate_figure(o4,
                       bottom = text_grob(" ", vjust=0.6,hjust=-1.5,size=13,face="bold"))+
  annotate("text", x=0.36,y=0.01,label="Accuracy Index: DR              Accuracy Index: DLT",size=3,color="#292929")+
  annotate("text", x=0.13,y=0.033,label=".",size=20,fontface="bold",color="#c6c88c")+
  annotate("text", x=0.15,y=0.033,label=".",size=20,fontface="bold",color="#E7E6E6")+
  annotate("text", x=0.36,y=0.033,label=".",size=20,fontface="bold",color="#736A35")+
  annotate("text", x=0.38,y=0.033,label=".",size=20,fontface="bold",color="#D0CECE")
mo4

pdf("drdlt_sec2_n100.pdf", width=7, height=10); mo4; dev.off()
