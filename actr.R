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
      print("Parameters of act-r verified")
    }
    return(actr_par)
}

## Base level activation computation.  Assumes a matrix "history" that has
## a column vector for each past retrieval moment, where the column vector
## identifies the winning item for each trial.

compute_base_levels <- function(moment,history,actr_par_retr) {

  #remove the failure from the history to not calculate an activation level for it
  history <- history[!history$chunk %in% "failure",]

  moments <- history$moment
  ## time since last retrieval for each retrieval (converted from milliseconds)
  tj <- (moment/1000) - (moments/1000);
  
  ## just pay attention to past, not future
  past <- history[ tj > 0,];
  
  ## the decay function
  tjd <- tj ^ (-actr_par_retr$d)
  decay <- tjd[tj > 0];
  
  ## compute base level activation for each chunk (for each trial)
  base_levels <- data.frame(chunk=unique(history$chunk));
  base_levels$levels <- NA
  for (c in unique(history$chunk)) {
    retrievals_boolean <- past$chunk == c;                  # boolean matrix
    activations <- retrievals_boolean * decay;
    b <- log(sum(activations,na.rm=TRUE));     # sum over all j retrievals
    b[is.infinite(b)] <- 0;                      # don't propagate through infinite values
    base_levels[base_levels$chunk ==c,]$levels <- b;
  }
  return(base_levels);
}



## Retrieval time computation: computes distribution of retrieval times for all
## items given a retrieval cue feature vector and a moment in time to do the
## retrieval.  Updates the history matrix with the winning items for this
## retrieval moment over all the trials.

retrieve <- function(cue_names,retrieval_cues,retrieval_moment,chunks,history,actr_par_retr) {

  num_features  <- length(cue_names)
  chunk_features <- chunks[,colnames(chunks) %in% cue_names ]
  creation_moment <- chunks$created
  num_chunks <- nrow(chunks)
  num_cues <- length(retrieval_cues[toupper(as.character(retrieval_cues))!="NULL"]);

    ## compute base level activations
    base_levels <- compute_base_levels(retrieval_moment,history,actr_par_retr)

    ## compute the match between chunks and retrieval cues (a boolean matrix)
    cues <- matrix(data=as.matrix(retrieval_cues),nrow=num_chunks, ncol=num_features,byrow=TRUE)
    is.nil <- chunk_features == "nil";    
    is_variable_cue <- cues == "VAR"; ##ASK FELIX

    match <- (chunk_features == cues);
    match_inc_var <- match | (is_variable_cue & !is.nil);

    ## checks which chunks exist at this moment
    exists <- matrix(creation_moment < retrieval_moment,nrow=num_chunks,ncol=num_features);

    ## checks which chunks match category cue
    chunk_category <- chunk_features[,cue_names=="cat"];
    cue_category <- cues[,cue_names=="cat"];
    matches_category <- chunk_category == cue_category

    ## compute fan for each feature: number of existing chunks matching each feature
    fan <- colSums(match_inc_var & exists) +  actr_par_retr$VAR.fan * is_variable_cue[1,];
    strength <- actr_par_retr$mas - log(fan)                       # fan equation


    ## compute source activation available for each cue (source spread over
    ## cues) and multiply by the fan (S * W in act-r equation).

    ## THIS IS NEW: We make VAR cues provide half the activation of other
    ## cues (because they only provide only half of the {feature,value} pair)
    cue.weights <- 1 - as.integer(is_variable_cue[1,])/2
    cue.weights <- cue.weights/sum(cue.weights[retrieval_cues!="NULL"])

    W <- actr_par_retr$G * cue.weights;    

#    W <- G/num_cues;  #ASK FELIX
    
    sw <- matrix(strength * W, nrow=num_chunks, ncol=num_features,byrow=TRUE);
    
    ## compute extra activation for each chunk; sum must ignore NA's because
    ## those correspond to features that are not retrieval cues.
    extra <- rowSums(match_inc_var * sw, na.rm=TRUE);    
    
    ## compute mismatch penalty
    is.retrieval.cue <- (cues != "NULL") & (cues != "VAR");

    if (actr_par_retr$var.mismatch.penalty) {
      mismatch <- (!match & is.retrieval.cue)  | (is_variable_cue & is.nil);
    } else {
      mismatch <- (!match & is.retrieval.cue);
    };

    ## mismatch <- (!match & is.retrieval.cue)  | (is_variable_cue & is.nil);
    penalty <- rowSums(mismatch * actr_par_retr$match.penalty);
    
    ## compute activation boost/penalty
    boost <- extra + penalty;
    

    ## add to base-level
    if (actr_par_retr$modulate.by.distinct) {
      ## compute how distinctive each chunk is (proportional to base-level activation)
      d.boost <- actr_par_retr$distinctiveness + base_levels$level;   
      activation <- base_levels$level + boost * d.boost;
    } else {
      activation <- base_levels$level + boost;
    };
    
    noise <- rlogis(num_chunks,0,actr_par_retr$ans) #
    noisy_activation <- activation + noise;

    ## make chunks that don't exist yet, or that don't match category cues,  have activation of -999
    exists <- creation_moment <  retrieval_moment

    nonexistence <- ifelse(exists,0,-Inf) #or NA?

    miscategory <- ifelse(matches_category,0,actr_par_retr$cat.penalty)

    final_activation  <- noisy_activation+nonexistence +miscategory

    final_activation <- c(final_activation,failure=actr_par_retr$rt)

    #the winner is the one with the highest activation if it's over 'rt'
    winner <- final_activation==max(final_activation)


    #to add a null word when no word gets to the retrieval treshold


    ## compute latency given the noisy activation, and mean and sd over all the
    ## monte carlo trials. Make non-existent chunks have a retrieval time of 9999.
    final_latency <- ifelse(winner,(actr_par_retr$F * exp(-final_activation))*1000,NA)

   


    

    result <- data.frame(chunk= names(final_latency),
      name= c(as.character(chunks$name),NA) ,
       final_latency=final_latency,
       final_activation=final_activation,
       winner=winner)


      update_history <- data.frame(chunk = result[result$winner==1,]$chunk,
        name = result[result$winner==1,]$name,
        moment=retrieval_moment+result[result$winner==1,]$final_latency )
  
    
    newhistory <- rbind(history,update_history )


return(list(output=result, newhistory=newhistory,finished_at=newhistory[nrow(newhistory),]$moment))
}


 


