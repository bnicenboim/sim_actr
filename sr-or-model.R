########################################################################################################
###
###  Simple Memory Retrieval Model of Agreement Phenomena (based on ACT-R theory)
### 
###  Bruno Nicenboim (bruno.nicenboim@uni-potsdam)
###
###  Reimplementation in R of ACT-R-in-R
###  orginal version https://github.com/felixengelmann/ACTR-in-R/blob/master/sr-or-model.r
### 
###  Original version edited by Felix Engelmann (felix.engelmann@uni-potsdam.de) 4. Oct. 2014:   
###  from an earlier script of 
###  Rick Lewis & William Badecker (rickl@umich.edu)
###  Version 3.0
###  2 Feb 2007
###
#########################################################################################################

rm(list=ls())
library(utils)
library(plyr)
library(lme4)
library(lme4)
library(ggplot2)
source("actr.r")
source("act-r-model.r")

RC <- list(

                  ## "the reporter who SENT the photographer to the editor hoped..."                  
                  SRC=list(
                       num_experimental_items = 10,
                       retrieval_schedule = "retrievals-subj-rel.txt", #or data frame or matrix
                       creation_schedule = "items-subj-rel.txt", #or data frame or matrix
                       procedural_duration = 100),

                  ## "the reporter who the photographer SENT to the editor hoped..."                  
                  ORC=list(
                       num_experimental_items = 10, 
                       retrieval_schedule = "retrievals-obj-rel.txt",
                       creation_schedule = "items-obj-rel.txt",
                       procedural_duration = 150)
                  )



#run_model will run nsims x num_experimental_items, while every row in actr_par is treated as a different subject

sim_RC <- sim_actr(RC,nsim=1)

summary(sim_RC,removeNaN=T)
summary(sim_RC,latencies=TRUE)

plot(sim_RC)



### Assuming individual differences in WMC:


actr_Gvar <- actr_default
actr_Gvar$G <- sort(rnorm(80,1,0.25))
actr_Gvar_subjs<- expand.grid(actr_Gvar)


sim_RC_G <- sim_actr(RC, actr_par=actr_Gvar_subjs,nsim=1)

summary(sim_RC_G,removeNaN=T)
summary(sim_RC_G,latencies=TRUE)


sim_data <- summary(sim_RC_G,latencies=TRUE,by_subj=T)

subjs <- actr_Gvar_subjs
subjs$subj <- seq_len(nrow(subjs))
sim_data <- merge(sim_data,subjs)

sim_data$G_group <- ifelse(scale(sim_data$G)>0,"high","low")

ddply(sim_data,.(G_group,retrieval_at,.id),summarize,mean(mean_latency))


sim_data$condition <- factor(sim_data$.id)
contrasts(sim_data$condition) <- contr.sum(2)

summary(m<-lmer(log(mean_latency)~condition *scale(G)+ (1|subj),data=sim_data,subset=retrieval_at=="VP2" ))

p<- ggplot(sim_data,aes(x=G,y=mean_latency,color=condition,linetype=condition, group=condition))+ facet_grid(. ~ retrieval_at) + geom_point() + geom_smooth()
print(p)

