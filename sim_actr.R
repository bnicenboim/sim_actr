library(plyr)

run_actr_each_condition <- function(data_list,actr_par=NULL,nsim=1) {
   
   if (is.null(actr_par)){
    actr_par <- actr_default_1subj
   }

  print("Starting to run model....")
  print(deparse(substitute(data_list)))
  print("****************")
  ## Read in the creation schedule:
  chunks <- create_schedule(data_list$creation_schedule)

  chunks$created <- as.numeric(as.character(chunks$created))
  
  #rownames word 1, word 2,....
  rownames(chunks) <- paste("word",1:nrow(chunks),sep="")

  #chunks
    #       name created cat mom_cat role number NOM ACC
    # word1  NP1       0  NP      IP spec   sing yes yes
    # word2  DP2     400  DP      CP spec    nil nil nil
    # word3  VP2     600  VP   I-bar comp   sing nil nil
    # word4  NP2     800  NP   I-bar spec   sing  no yes
    # word5  NP3    1400  NP      VP comp   sing  no yes
    # word6  VP1    1800  VP      IP comp   sing nil nil

  retrievals <- create_schedule(data_list$retrieval_schedule)
    #each retrieval is names retr1...
  retrievals$retr <- paste("retr",1:nrow(retrievals),sep="")
#   retrievals
#   moment cat mom_cat role number  NOM  ACC  retr
# 1    DP2  NP      IP spec   sing  yes NULL retr1
# 2    VP2  DP      CP spec   NULL NULL NULL retr2
# 3    VP1  NP      IP spec   sing  yes NULL retr3

  cue_names <- colnames(chunks[-(1:2)])
  cue_names_retr <- colnames(retrievals[ -c(1,ncol(retrievals)) ])
  
  #check if retrieval and creation have the same cues
  #(catch if there's a typo)
  if(!all(cue_names %in% cue_names_retr) | !all(cue_names_retr %in% cue_names)  ) {
    warning("Cues from the retrieval and creation schedule are not identical")
  }



  #verify the actr params
  actr_par <- verify_actr_par(actr_par)

  #simulation of trials
  #every row in actr_par can be treated as a different subject
  nsubj <- nrow(actr_par)
  nretrievals <- nrow(retrievals)
  unique_trials <-   data.frame(exp= rep(1:nsim,each=nsubj*data_list$num_experimental_items ),
    subj=rep(1:nsubj, each=data_list$num_experimental_items,times=nsim),
    item=rep(1:data_list$num_experimental_items,times=nsim*nsubj)
    )
  unique_trials_params <- cbind(unique_trials,repeat_row(actr_par,each=data_list$num_experimental_items,times=nsim))

  rep_trials <- (repeat_row(unique_trials_params,each=nretrievals))
  rownames(rep_trials)<- NULL
  rep_retr <- repeat_row(retrievals,times=nsubj*data_list$num_experimental_items*nsim)
  rownames(rep_retr) <- NULL
 
  #all the retrievals schedule in all the trials, moment is at which creation the retrieval is triggered:
  sim_retrievals <- cbind(rep_trials,rep_retr)  

 
  nunique_trials<- nrow(unique_trials)

  trials_history <- repeat_row(cbind(wordn=rownames(chunks),chunks),times=nunique_trials)
  #I add a new level failure for when no chunk is retrieved:
  trials_history$name <- factor(trials_history$name, levels=c(levels(chunks$name),"failure"))
  
  trials_h <- repeat_row(unique_trials,each=nrow(chunks))
  rownames(trials_h) <- NULL
  history <- cbind(trials_h,trials_history)


  available_moment <- cbind(unique_trials,available_at=-1) #before any retrieval
  
  repeated_par <- repeat_row(actr_par, each=data_list$num_experimental_items)  
  rownames(repeated_par) <- NULL
  retrieval_moment <- cbind(unique_trials,repeated_par)
  retrieval_moment$moment <- NA
  res <- data.frame()
  
  for(r in unique(sim_retrievals$retr)){
    #this checks if the retrieval can be done as schedule (a[r,]$moment) or delayed to the next available moment (available_moment) #that is, after the previous retrieval was completed
    scheduled_moment  <- chunks[chunks$name %in% sim_retrievals[sim_retrievals$retr %in% r,]$moment,]$created

    retrieval_moment$moment <- pmax(available_moment$available_at, scheduled_moment )

    retrieval_cues <- sim_retrievals[sim_retrievals$retr %in% r,colnames(sim_retrievals) %in% cue_names]

    #all the retrieval cues should be the same, since these are different trials of the same retrieval
    if(nrow(unique(retrieval_cues))!=1){
      stop("something's wrong!!")
    }
    retrieval_cues <- unique(retrieval_cues)

    history <- history[order(history$exp,history$subj,history$item,history$created),]

    #retrieval using act-r formulas here:
    ret <- retrieve(cue_names,retrieval_cues, retrieval_moment,history)

    history <- ret$newhistory
    ret$output$retr <- r

    ret$output$retrieval_at <- as.character(sim_retrievals[sim_retrievals$retr %in% r,]$moment)
    available_moment$available_at <- ret$finished_at

    #ASK FELIX IF I SHOULD LEAVE TIME FOR THE RULE
    res <- rbind(res, ret$output)
    }
    
    
  results<-res    

  #old ddply:
  #latencies <- ddply(subset(results,winner==1),.(exp,subj,item,retrieval_at),summarize,latency=sum(final_latency))

  latencies <- summarise(group_by(subset(results,winner==1),exp,subj,item,retrieval_at),latency=sum(final_latency))



  complete_results =list(results=results,latencies=latencies,history=history)
  return(complete_results)
}



sim_actr <- function(conditions_list,actr_par=NULL,nsim=1){
  if (is.null(actr_par)){
    print("Using default act-r parameters. See actr_default")
    actr_par <- actr_default_1subj
  }

results <- llply(conditions_list, run_actr_each_condition, actr_par=actr_par,nsim=nsim)

output <- list()
output$results <- results
output$actr_par <- actr_par
output$nsim <- nsim

#I don't know what's better
#class(output) <- append(class(output),"actr")
attr(output, "class") <- "actr"

  return(output)  
}



repeat_row <- function(data,...){
  data[ rep(seq_len(nrow(data)), ...),]  
}


actr_default <- list(
## Latency factor
 F = c(0.14),
## Extra category penalty
cat.penalty = c(-Inf),
## Total source activation
G = c(1.0),  
## Activation noise parameter for logistic distribution
ans  = c(0.15), #lewis vasishth?
## Fan parameter
mas = c(1.5), #lewis vasishth?
## Base level decay parameter
d = c(0.5), #standard
## Match penalty
# match.penalty <- c(0)
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
rt = c(-1.5) 
)

actr_default_1subj <- as.data.frame(actr_default)



create_schedule <- function(unformated_schedule){

  if(is.data.frame(unformated_schedule)){
    schedule <- unformated_schedule
  }  else if(is.matrix(unformated_schedule)){
    schedule <- data.frame(unformated_schedule)
  } else if(file.exists(unformated_schedule)){
    schedule <- read.delim(file=unformated_schedule,header=FALSE,colClasses="character", row.names=1)
  } else {
    stop(paste("Error in ",names(args)[1]," file or data.frame"))
  }

  #This is to flip the data frame in case moment/name cat etc are not the colnames
  if(!"cat" %in% colnames(schedule)){
    schedule <- t(schedule)
    rownames(schedule) <- NULL
    schedule<-data.frame(schedule)
  }
  
  #validation of the row names
  colnames(schedule) <- make.names(colnames(schedule))
  return(schedule)
  
}



summary.actr <- function(output,latencies=FALSE,by_subj=FALSE,pars=NULL,removeNaN=FALSE,abbreviate=TRUE){

  #set the right groups for the dplyr later:
  pars_subj <- NULL
  all_par <- NULL
  if(by_subj){
    if(!is.null(pars)){
      all_par <- output$actr_par[colnames(output$actr_par) %in% pars]
      all_par$subj <- seq_len(nrow(all_par))
    }

    pars_subj <- c(pars,"subj")
  }    
  

  summ <-  ldply(output$results, function(l){
  
    df <- data.frame()

    if(latencies){
      by <- c(pars_subj,"retrieval_at")
      # Convert character vector to list of symbols
      dots <- lapply(by, as.symbol)

      if(!is.null(all_par)){
        latencies <- left_join(l$latencies,all_par)
        grouped <- group_by_(latencies,.dots=dots)
      } else {
        grouped <- group_by_(l$latencies,.dots=dots)
      }

      df<- summarise(grouped,Latency = mean(latency), SE= sd(latency)/sqrt(length(latency)))

    } else if(!latencies){
      by <- c(pars_subj,"retr","wordn","name")
      dots <- lapply(by, as.symbol)

      if(!is.null(all_par)){
        results <- left_join(l$results,all_par)
        grouped <- group_by_(results,.dots=dots)
  
      } else {
        grouped <- group_by_(l$results,.dots=dots)
      }

      df <- summarise(grouped,Retr_P=mean(winner),Latency=mean(final_latency,na.rm=T),mean_activation=mean(final_activation))

      if(abbreviate){
        df$name <- abbreviate(iconv(df$name, to='ASCII//TRANSLIT'),12,"both")
      }
    }
    return(df)
  }) #end of ldply


  if(removeNaN){
    summ <- summ[!is.nan(summ$Latency),]
  }
  
  if(abbreviate){
    summ$.id <- abbreviate(iconv(summ$.id, to='ASCII//TRANSLIT'),12,"both")
  }

  nitem<-  llply(output$results, function(l){
    return(max(l$results$item))
  })


  lines <- c("# Summary of simulations",
    paste("# Number of simulations: ",output$nsim),
    paste(paste("# Number of items ",names(nitem)),nitem),
    paste("# Number of simulated subjects: ",nrow(output$actr_par),"\n"))

  cat(paste(lines,collapse="\n"))

  return(summ)
}


plot.actr <- function(output,pars=NULL){

  by_subj <- ifelse(is.null(pars),FALSE,TRUE)
  summary <- summary(output,latencies=TRUE,pars=pars[1],by_subj=by_subj) #only the first par

  if(is.null(pars)){
    plot <- ggplot(summary, aes(x=.id, y=Latency, fill=.id))+ 
            facet_grid(. ~ retrieval_at) + 
            geom_bar(position=position_dodge(), stat="identity") +
            geom_errorbar(aes(ymin=Latency-2*SE, ymax=Latency+2*SE),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))+
            xlab("Experimental Condition") +
      ylab("Latency") +
      scale_fill_hue(name="Experimental Condition") +
      ggtitle("Summary") +
      #scale_y_continuous(breaks=seq(0:summary$)) +
      theme_bw()
  } else {
    plot <- ggplot(summary, aes_string(x=pars[1], y="Latency", color=".id",linetype=".id",shape=".id"))+ 
            facet_grid(. ~ retrieval_at) + facet_grid(. ~ retrieval_at) + geom_point() + geom_smooth()+

            xlab(paste("Value of act-r param",pars[1])) +
      ylab("Latency")  +
      ggtitle("Summary") +
      #scale_y_continuous(breaks=seq(0:summary$)) +
      theme_bw()
  }



  return(plot)
}
