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
#library(plyr)
library(dplyr)
library(lme4)
library(lme4)
library(ggplot2)
source("actr.R")
source("sim_actr.R")

read_data_frame <- function(text){
  #Convenience function to read tables in Hadley Wickham's data_frame
  as_data_frame(read.table(text=text,header=T,stringsAsFactors=F))
}

#the reporter WHO SENT the photographer to the editor HOPED...
SRC <- read_data_frame("
name             created cat mom_cat role number case
the_reporter     0       NP  IP      spec sing   nom
WHO              450     DP  CP      spec sing    nom
SENT             650     VP  I-bar   comp sing   nil
the_photographer 800     NP  I-bar   spec sing   acc
to_the_editor    1400    NP  VP      comp sing   acc
HOPED            1850    VP  IP      comp sing   nil
")


retr_SRC<- read_data_frame("
moment           cat mom_cat role number case
WHO              NP  IP      spec sing   nom
SENT             DP  CP      spec sing   nom
the_photographer VP  I-bar   comp sing   nil
to_the_editor    VP  I-bar   comp sing   nil
HOPED            NP  IP      spec sing   nom
")


#"the reporter WHO the photographer SENT to the editor HOPED..."
ORC <- read_data_frame("
name             created cat mom_cat role number case
the_reporter     0       NP  IP      spec sing   nom
WHO              450     DP  CP      spec sing    acc
the_photographer 600     NP  I-bar   spec sing   nom
SENT             1050    VP  I-bar   comp sing   nil
to_the_editor    1450    NP  VP      comp sing   acc
HOPED            1900    VP  IP      comp sing   nil
")



retr_ORC<-read_data_frame("
moment           cat mom_cat role number case
WHO              NP  IP      spec sing   nom
SENT             DP  CP      spec sing   acc
SENT             NP  I-bar   spec sing   nom
to_the_editor    VP  I-bar   comp sing   nil
HOPED            NP  IP      spec sing   nom
")



RC <- list(

                  ## "the reporter WHO SENT the photographer to the editor HOPED..."                  
                  SRC=list(
                       num_experimental_items = 10,
                       creation_schedule = SRC, 
                       retrieval_schedule = retr_SRC, 
                       procedural_duration = 50), 

                  ## "the reporter WHO the photographer SENT to the editor HOPED..."                  
                  ORC=list(
                       num_experimental_items = 10, 
                       creation_schedule = ORC,
                       retrieval_schedule = retr_ORC,
                       procedural_duration = 150) 
                  )



#run_model will run nsims x num_experimental_items, while every row in actr_par is treated as a different subject


sim_RC <- sim_actr(RC,nsim=100)

print(summary(sim_RC))
print(summary(sim_RC),removeNaN=T)


plot(sim_RC)





### Assuming individual differences in WMC:


actr_Gvar <- actr_defaults
actr_Gvar$G <- sort(rnorm(80,1,0.25))
actr_Gvar_subjs<- as_data_frame(expand.grid(actr_Gvar))

#assumes that everyline of actr_Gvar_subjs is comming from a different subject
head(actr_Gvar_subjs)

sim_RC_G <- sim_actr(RC, actr_par=actr_Gvar_subjs,nsim=1)

print(summary(sim_RC_G),minlength=8)

plot(sim_RC_G,pars="G")


#other analysis:
sim_data <- by_subj(sim_RC_G,pars="G")

sim_data$G_group <- ifelse(scale(sim_data$G)>0,"high","low")

summarise(group_by_(sim_data,
  .dots=c("G_group","retrieval_at",".id"))
,mean(Latency))


sim_data$condition <- factor(sim_data$.id)
contrasts(sim_data$condition) <- contr.sum(2)

summary(m<-lmer(log(Latency)~condition *scale(G)+ (1|subj),data=sim_data,subset=retrieval_at=="SENT" ))

