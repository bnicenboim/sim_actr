verify_actr_par <- function(actr_par){
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

    model_par <- names(actr_par) %in% names(actr_defaults)
    if(!all(model_par)){
      irrelevant_par <- colnames(actr_par)[which(!model_par)]
      warning(paste("Unknown parameters in the model:",paste(irrelevant_par,collapse="; " )))
    }
    
    found_par <- names(actr_defaults) %in% names(actr_par)

       if(!all(found_par)){
        missing_par <- names(actr_defaults)[which(!found_par)]
        warning(paste("Missing parameters in the model (using default values):",paste(missing_par,collapse="; " )))

      actr_par[eval(missing_par)] <- actr_defaults[eval(missing_par)]
    }

    if(all(model_par) & all(found_par)){
      #print("Parameters of act-r verified")
    }
    return(actr_par)
}

## Base level activation computation.  Assumes a matrix "history" that has
## a column vector for each past retrieval moment, where the column vector
## identifies the winning item for each trial.

compute_base_levels <- function(retrieval_moment,history) {

  #remove the failure from the history to not calculate an activation level for it
  chistory <- history[!history$name %in% "failure",]
  chunk_names <- unique(chistory$name)
  nchunks <- length(chunk_names)
  ## time since last retrieval for each retrieval (converted from milliseconds)
  
  #chistory <- merge(chistory,retrieval_moment[c("exp","subj","item","moment","d")])

  #new dplyr
  chistory <- inner_join(chistory,retrieval_moment[c("exp","subj","item","moment","d")],by=c("exp", "subj", "item"))


  chistory$tj <- (chistory$moment - chistory$created)/1000
  
  
  ## just pay attention to past, not future
  chistory[ chistory$tj <= 0,]$tj <- NA

  
  ## the decay function
  chistory$activation <- chistory$tj ^ -(chistory$d)

  #old ddply:  
  #  base_levels<- ddply(chistory,.(exp,subj,item,wordn,name),summarize, activation=log(sum(activation,na.rm=TRUE)))

  base_levels<- summarise(group_by(chistory, exp,subj,item,wordn,name),
    activation=log(sum(activation,na.rm=TRUE)))

  
  base_levels$activation <- ifelse(is.infinite(base_levels$activation),0,base_levels$activation)

  return(base_levels);
}



## Retrieval time computation: computes distribution of retrieval times for all
## items given a retrieval cue feature vector and a moment in time to do the
## retrieval.  Updates the history matrix with the winning items for this
## retrieval moment over all the trials.

retrieve <- function(cue_names,retrieval_cues,retrieval_moment,history) {

    ## compute base level activations
  base_levels <- compute_base_levels(retrieval_moment,history)
  
  nchunks <- length(unique(base_levels$name))

  cues <- unique(history[colnames(history) %in% c("name",cue_names)])
  #dplyr
  #base_levels_cues <- merge(base_levels,cues,by="name",all.x=TRUE,all.y=FALSE,sort=FALSE)
  base_levels_cues <- left_join(base_levels,cues,by="name")
  base_levels_cues <- base_levels_cues[order(base_levels_cues$exp,base_levels_cues$subj,base_levels_cues$item,base_levels_cues$wordn),] #previous order

  trial_info <-base_levels_cues[c("exp","subj","item")]
  unique_trials <- unique(trial_info)

  base_levels_cues_m <- as.matrix(base_levels_cues[colnames(base_levels_cues) %in% cue_names])
  

  retrieval_cues_m <- as.matrix(repeat_row(retrieval_cues,each=nrow(base_levels_cues)))

  retrieval_moment_d <- repeat_row(retrieval_moment,each=nchunks)

  match <- base_levels_cues_m == retrieval_cues_m

  #nil are irrelevant features?
  is_nil <- base_levels_cues_m == "nil"

  #VAR is cue that's going to match any feature that exists (not nil)
  is_variable_cue <- retrieval_cues_m =="VAR" 
  match_inc_var <- match | (is_variable_cue & !is_nil)

   ## checks which chunks exist at this moment
  if(0){
   earliest_history <-  summarise(group_by(subset(history,!name %in% "failure"), exp,subj,item,wordn,name),earliest_moment= min(created))
  

    exists <- earliest_history$earliest_moment < retrieval_moment_d$moment
  }
  # SEEE IF IT WORKS!!!!!!!!
  exists <- !is.infinite(base_levels$activation) & !is.na(base_levels$activation) &  !is.nan(base_levels$activation)

  matches_category <- base_levels_cues_m[,"cat"] == retrieval_cues_m[,"cat"]

  match_inc_var_exist <- match_inc_var & exists


  match_inc_var_exist_trial_info<- cbind(trial_info,match_inc_var_exist)
  #old ddply
  #fan <- ddply(match_inc_var_exist_trial_info,.(exp,subj,item),colwise(sum,colnames(retrieval_cues)))

  #we calculate the fan for each cue ignoring the features of other cues (values in columns). In real act-r, the fan is computed over all slots in all chunks 
  fan <- summarise_each(group_by(match_inc_var_exist_trial_info,exp,subj,item),funs(sum))
  print(fan)

  #this is some Rick Lewis thing:
  fan_update <- repeat_row(retrieval_cues=="VAR",each=nrow(unique_trials)) * retrieval_moment$VAR.fan

  fan[colnames(fan)%in% cue_names]<- fan[colnames(fan)%in% cue_names] +fan_update

  #maximum associated strenght(maximum of spreading activation that can go from a cue to a chunk) - log of fan => the bigger the fan, => less strength
  strength <- retrieval_moment$mas - log(fan[colnames(fan)%in% cue_names])

  strength[strength==Inf] <- NA


  ## compute source activation available for each cue (source spread over
  ## cues) and multiply by the fan (S * W in act-r equation).
  ## THIS IS NEW: We make VAR cues provide half the activation of other
  ## cues (because they only provide only half of the {feature,value} pair)
  cue_weights_temp <- 1 - (is_variable_cue)/2
 
  #ASK FELIX , WHY DO I DIVIDE BY SUM?????
  #cue_weights <- cue_weights/sum(cue_weights[retrieval_cues!="NULL"])
  #SHOULD IT BE LIKE THIS:
 # cue_weights_temp_trial <- cbind(trial_info,cue_weights_temp)
 no_nulls <- colnames(retrieval_cues[,retrieval_cues!="NULL"])

 # div_cue_weights_temp <- ddply(cue_weights_temp_trial,.(exp,subj,item),colwise(sum,no_nulls))
 # div_cue_weights_temp$div <- rowSums(div_cue_weights_temp[no_nulls])
  
  # cue_weights <- cue_weights_temp/rep(   div_cue_weights_temp$div,each=nchunks)
  cue_weights <- cue_weights_temp/rowSums(cue_weights_temp[,no_nulls])

  W <- retrieval_moment_d$G  * cue_weights #= G/number of relevant features (if VAR isn't set)

  strength_m <- repeat_row(strength,each=nchunks)
  
  #final amount of spreading activation from cue to features:
  sw <- strength_m  * W

  #where it doesn't match there's no spreading activation:
  extra <- rowSums(match_inc_var * sw, na.rm=TRUE)    

  is_retrieval_cue <- (retrieval_cues_m != "NULL") & (retrieval_cues_m != "VAR");

  #check the var thing again
  #TRUE is the cues which is relevant and doesn't match
  mismatch <- !match & is_retrieval_cue | 
              (is_variable_cue & is_nil) & retrieval_moment_d$var.mismatch.penalty #this line is relevant if var.mismatch.penalty is T




  penalty <- rowSums(mismatch *retrieval_moment_d$match.penalty)
    ## compute activation boost/penalty
  
  boost <- extra + penalty;
    
  ## compute how distinctive each chunk is (proportional to base-level activation)
  d.boost <- ifelse(retrieval_moment_d$modulate.by.distinct, retrieval_moment_d$distinctiveness + base_levels$activation ,1)

  activation <- base_levels$activation + boost * d.boost #d.boost=1 unless it's modulated by distinct
    
  noise <- rlogis(nrow(base_levels),0,retrieval_moment_d$ans) #
  noisy_activation <- activation + noise;

    ## make chunks that don't exist yet, or that don't match category cues,  have activation of -Inf
  nonexistence <- ifelse(exists,0,-Inf) #or NA?

  miscategory <- ifelse(matches_category,0,retrieval_moment_d$cat.penalty)

  final_activation  <- noisy_activation+nonexistence +miscategory

  base_levels$final_activation <- final_activation
  base_levels$F <- retrieval_moment_d$F
  #omission errors:
  omissions<-retrieval_moment[c("exp","subj","item")]
  omissions$wordn <-NA  
  omissions$name <- "failure"
  omissions$activation <-NA
  omissions$final_activation <- retrieval_moment$rt
  omissions$F <- retrieval_moment$F #I'll need it for latencies later

  base_levels_output <- rbind.fill(base_levels,omissions)


    #the winner is the one with the highest activation if it's over 'rt'
   
   #old ddply
    # winner<- ddply(base_levels_output,.(exp,subj,item), function(df){
    #  retrieved = df[which(max(df$final_activation)==df$final_activation),]$name
    #  return(data.frame(retrieved=retrieved))
    #   })
  winner <- filter(group_by(base_levels_output,exp,subj,item),  max(final_activation) ==final_activation)
    winner$retrieved <- winner$name
    winner<-winner[c("exp","subj","item","retrieved")]

  #dplyr
  #base_levels_output <- merge(base_levels_output,winner,all.x=T)
  base_levels_output <- left_join(base_levels_output,winner,by=c("exp", "subj", "item"))
  base_levels_output$winner <- base_levels_output$retrieved == base_levels_output$name

    

  ## compute latency given the noisy activation

  base_levels_output$final_latency <- ifelse(base_levels_output$winner,(base_levels_output$F * exp(-base_levels_output$final_activation))*1000,NA)

   


  result <- base_levels_output[c("exp","subj","item","wordn","name","final_activation","winner","final_latency")]

 #dplyr
 # update_winner<-merge(result[result$winner,],cues,by=c("name"),all.x=T)
  update_winner<-left_join(result[result$winner,],cues,by=c("name"))
  
  update_winner$created <- retrieval_moment$moment +update_winner$final_latency


  update_history <- rbind.fill(history,update_winner )
  newhistory <- update_history[,colnames(history)]

return(list(output=result, newhistory=newhistory,finished_at=update_winner$created))
}


 
