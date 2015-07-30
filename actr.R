#defaults:
  actr_defaults <- list(
    ## Latency factor
     F = c(0.14),
    ## Extra category penalty
    cat.penalty = c(-Inf),
    ## Total source activation
    G = c(1.0),  
    ## Activation noise parameter for logistic distribution
    ans  = c(0.15), 
    ## Fan parameter
    mas = c(1.5), #lewis vasishth?
    ## Base level decay parameter
    d = c(0.5), #standard
    ## Match penalty
    match.penalty = c(-1.5),
    ## VAR mismatch penalty
    var.mismatch.penalty = c(FALSE),
    ## Additional fan associated with VAR retrieval cues
    VAR.fan = c(0),
    ## Distinctiveness parameters.  Note that the quantitative parameter has no
    ## effect when modulate.by.distinct is FALSE. So this is some  wasted
    ## effort in the simple code below.
    modulate.by.distinct = c(FALSE),
    distinctiveness = c(0),
    #retrieval treshold
    rt = c(-1.5) #-1.5
    )
actr_par_names <- names(actr_defaults)

verify_actr_par <- function(actr_par){
  

    model_par <- names(actr_par) %in% actr_par_names
    if(!all(model_par)){
      irrelevant_par <- colnames(actr_par)[which(!model_par)]
      warning(paste("Unknown parameters in the model:",paste(irrelevant_par,collapse="; " )))
    }
    
    found_par <- actr_par_names %in% names(actr_par)

       if(!all(found_par)){
        missing_par <- actr_par_names[which(!found_par)]
        warning(paste("Missing parameters in the model (using default values):",paste(missing_par,collapse="; " )))

      actr_par[eval(missing_par)] <- actr_defaults[eval(missing_par)]
    }

    if(all(model_par) & all(found_par)){
      #print("Parameters of act-r verified")
    }
    return(actr_par)
}

## Base level activation computation.  Assumes a matrix "history_till_now" that has
## a column  for each past retrieval moment

compute_base_levels <- function(history_till_now,cue_names) {

  #remove the failure from the history to not calculate an activation level for it
  historic_chunks <- filter(history_till_now,name != "<failure>")
  
  ## time since last retrieval for each retrieval (converted from milliseconds)
  historic_chunks$tj <- (historic_chunks$now - historic_chunks$accessed)/1000
  
  ## the decay function
  historic_chunks$activation <- historic_chunks$tj ^ -(historic_chunks$d)
  
  # col_base_levels <- colnames(historic_chunks)[!colnames(historic_chunks) %in% c("moment","tj","activation")]

  col_base_levels <- select_vars_(names(historic_chunks),c("exp","subj","item","name","chunkn","now",cue_names,actr_par_names))
  
  base_levels<-summarise(  group_by_(historic_chunks,   .dots=col_base_levels),
    activation=log(sum(activation)))

unique(base_levels[c("exp","subj","item")])
  
   return(base_levels);
}



## Retrieval time computation: computes distribution of retrieval times for all
## items given a retrieval cue feature vector and a moment in time to do the
## retrieval.  Updates the history matrix with the winning items for this
## retrieval moment over all the trials.

retrieve <- function(cue_names,retrieval_cues,retrieval_moment,history_till_now) {
  
  ## compute base level activations
  base_levels <- compute_base_levels(history_till_now,cue_names)
  base_levels <- group_by(ungroup(base_levels),exp,subj,item)
  #arrange(select(base_levels,activation,now, everything()),activation,exp,subj,item)
  
  #(pseudo) matrix of cues for the chunks with base levels
  base_levels_cues_m <- select_(ungroup(base_levels),.dots=cue_names)  

  #(pseudo) matrix of retrieval cues
  retrieval_cues_m <- repeat_row(retrieval_cues,each=nrow(base_levels))

  match <- base_levels[cue_names]==retrieval_cues_m

  #nil are irrelevant features?
  is_nil <- base_levels[cue_names] == "nil"

  #VAR is cue that's going to match any feature that exists (not nil)
  is_variable_cue <- retrieval_cues_m =="VAR" 
  match_inc_var <- match | (is_variable_cue & !is_nil)

  matches_category <- base_levels_cues_m[,"cat"] == retrieval_cues_m[,"cat"]

  
  match_base_levels<- base_levels
  match_base_levels[cue_names] <- match_inc_var
  #old ddply
  #we calculate the fan for each cue ignoring the features of other cues (values in columns). In real act-r, the fan is computed over all slots in all chunks 
  fan <- summarise_each_(group_by(match_base_levels,exp,subj,item,VAR.fan,mas),funs(sum),cue_names) #so the there will be one row for each trial (subj x exp x item ) compare with base_levels (subj x exp x item x created chunks)

  #this is some Rick Lewis thing:
  fan_update <- repeat_row(retrieval_cues=="VAR",each=nrow(fan)) * fan$VAR.fan

  fan[cue_names]<- fan[colnames(fan)%in% cue_names] +fan_update

  #maximum associated strenght(maximum of spreading activation that can go from a cue to a chunk) - log of fan => the bigger the fan, => less strength
  strength <- fan #it has the same structure as fan:
  strength[cue_names] <- fan$mas - log(fan[cue_names])

  strength[strength==Inf] <- NA


  ## compute source activation available for each cue (source spread over
  ## cues) and multiply by the fan (S * W in act-r equation).
  ## THIS IS NEW: We make VAR cues provide half the activation of other
  ## cues (because they only provide only half of the {feature,value} pair)
  
  #TODO: a cleaner version with dplyr
  cue_weights_temp <- 1 - (is_variable_cue)/2 #nrow as baselevel
  no_nulls <- colnames(retrieval_cues[,retrieval_cues!="NULL"])
  cue_weights <- cue_weights_temp/rowSums(cue_weights_temp[,no_nulls])
  
  W <- base_levels
  W[cue_names] <- base_levels$G  * cue_weights #= G/number of relevant features (if VAR isn't set)

  
  #increase the size of strength_m data frame to the base_levels size:
  strength_m <- left_join(select(W,exp,subj,item),strength,by = c("exp", "subj", "item"))
  
  if(nrow(strength_m)!=nrow(base_levels)) {
    stop("Something's wrong in the base_levels data_frame")
  }

  #final amount of spreading activation from cue to features:
  sw <- strength_m[cue_names]  * W[cue_names]

  #where it doesn't match there's no spreading activation:
  extra <- rowSums(match_inc_var * sw, na.rm=TRUE)    

  is_retrieval_cue <- (retrieval_cues_m != "NULL") & (retrieval_cues_m != "VAR");

  #check the var thing again
  #TRUE is the cues which is relevant and doesn't match
  mismatch <- !match & is_retrieval_cue | 
              (is_variable_cue & is_nil) & base_levels$var.mismatch.penalty #this line is relevant if var.mismatch.penalty is T




  penalty <- rowSums(mismatch *base_levels$match.penalty)
    ## compute activation boost/penalty
  
  boost <- extra + penalty;
    
  ## compute how distinctive each chunk is (proportional to base-level activation)

  d.boost <- ifelse(base_levels$modulate.by.distinct, base_levels$distinctiveness + base_levels$activation ,1)

  activation <- base_levels$activation + boost * d.boost #d.boost=1 unless it's modulated by distinct
    
  noise <- rlogis(nrow(base_levels),0,base_levels$ans) #
  noisy_activation <- activation + noise;

  miscategory <- ifelse(matches_category,0,base_levels$cat.penalty)

  final_activation  <- noisy_activation +miscategory
  
  #this is irrelevant since what counts now is the final_activation
  base_levels <- select(base_levels,-activation)

  base_levels$final_activation <- final_activation
  #omission errors:
  omissions<-unique(base_levels[c("exp","subj","item",actr_par_names,"now")])
  omissions$chunkn <-"<failure>" 
  omissions$name <- "<failure>"
#  omissions$activation <-NA
  

  #the final activation in the case of omission is based on rt 
  omissions$final_activation <- omissions$rt

  base_levels_output <- bind_rows(base_levels,omissions)


  #the winner is the one with the highest activation    

  
  winner <- filter(group_by(base_levels_output,exp,subj,item),  max(final_activation) == final_activation)

  winner$retrieved <- winner$name
  
  winner<-winner[c("exp","subj","item","retrieved")]

  #dplyr
  #base_levels_output <- merge(base_levels_output,winner,all.x=T)
  base_levels_output <- left_join(base_levels_output,winner,by=c("exp", "subj", "item"))
  base_levels_output$winner <- base_levels_output$retrieved == base_levels_output$name

    

  ## compute latency given the noisy activation

  base_levels_output$final_latency <- ifelse(base_levels_output$winner,(base_levels_output$F * exp(-base_levels_output$final_activation))*1000,NA)

   
  base_levels_output <- select(base_levels_output,exp,subj,item,name,chunkn,winner,final_activation,final_latency,retrieved,everything())
  



  winners<- base_levels_output[base_levels_output$winner==1,]
  
  winners$accessed <- winners$now +winners$final_latency

  #irrelevant columns:
  #retrieved is irrelevant because it will be the same as name ()
  winners <- select(winners,-retrieved,-winner,-final_latency)
  
  newhistory_till_now <- bind_rows(history_till_now,winners)

return(list(actr_sim=base_levels_output, newhistory_till_now=newhistory_till_now,winners=winners))
}


 
