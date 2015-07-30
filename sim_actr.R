
#detach(package:plyr)
#it can create problems with dplyr
library(tidyr)
library(dplyr)

run_actr_each_condition <- function(data_list,actr_par=NULL,nsim=1,debug=F) {
   #For TEST only with
   #data_list <- RC$ORC
   #actr_par <-  repeat_row(actr_default_1subj,each=2)
   #actr_par[1,] <- actr_par[1,] +.5
   #nsim <- 3

   if (is.null(actr_par)){
    actr_par <- actr_default_1subj
   }
   actr_par <- as_data_frame(actr_par)

  print("Starting to run model....")
  print(deparse(substitute(data_list)))
  print("****************")
  ## Read in the creation schedule:
  chunks <- data_list$creation_schedule
  
 

  retrievals <- data_list$retrieval_schedule
    #each retrieval is names retr1...

  #all this to take all the column that are left from created
  cue_names <- colnames(chunks[(which(colnames(chunks)=="created")+1):ncol(chunks)])
  #everything after moment is a retrieval cue
  retr_cue_names <- colnames(select(retrievals,-moment ))
  
  #check if retrieval and creation have the same cues
  #(catch if there's a typo)
  if(!all(cue_names %in% retr_cue_names) | !all(retr_cue_names %in% cue_names)  ) {
    warning("Cues from the retrieval and creation schedule are not identical")
  }

  sanity_check(chunks,retrievals)

  nchunks <- nrow(chunks)
  #rownames chunk 1, chunk 2,....
  chunks$chunkn <- paste("chunk",1:nchunks,sep="")
  retrievals$retr <- paste("retr",1:nrow(retrievals),sep="")

 


  #verify the actr params
  actr_par <- verify_actr_par(actr_par)

  #simulation of trials
  #every row in actr_par can be treated as a different subject
  nsubj <- nrow(actr_par)
  nretrievals <- nrow(retrievals)
  nunique_trials<-nsim*nsubj*data_list$num_experimental_items

  #Build the trials data frame
  unique_trials <-   data_frame(
    exp= rep(1:nsim,each=nsubj*data_list$num_experimental_items),
    subj=rep(1:nsubj, each=data_list$num_experimental_items,times=nsim),
    item=rep(1:data_list$num_experimental_items,times=nsim*nsubj)
    )
  
  #subj_indiv <- repeat_row(actr_par,each=data_list$num_experimental_items,times=nsim)
  subj_indiv <-actr_par
  subj_indiv$subj <- 1:nrow(subj_indiv)

  rep_trials <- repeat_row(unique_trials,each=nretrievals)

  #each retrieval will be calculated for each experiment(nsim) x subj x item
  rep_retr <- repeat_row(retrievals,times=nunique_trials)
 
  #all the retrievals schedule in all the trials, moment is at which creation the retrieval is triggered:
  #it seems to be a bug that I need to force it data_frame again
  retrievals_table <- tbl_df(bind_cols(rep_trials,rep_retr))



  rep_history <- repeat_row(chunks,times=nunique_trials)
  
  
  unique_trials_pars <- full_join(unique_trials,subj_indiv,by="subj")
  trials_history <- repeat_row(unique_trials_pars,each=nchunks)
  
  history <- as_data_frame(bind_cols(rep_history,trials_history))
  history$name <- factor(history$name)
  #I add a new level failure for when no chunk is retrieved:
  history$name <- factor(history$name, levels=c(levels(history$name),"<failure>"))
  
  history<- rename(history, accessed = created)
  #the column accessed will contain the creation time at the beggining but then the last time a chunk was retrieved

  


  available_moment <- unique_trials
  available_moment$available_at <- -1  #before any retrieval
  
  retrieval_moment <- unique_trials
  retrieval_moment$now <- NA

  
  result_list = list()
  

  for(r in 1:nrow(retrievals)){
    #this checks if the retrieval can be done as schedule (retrieval_moment$now) or delayed to the next available moment (available_moment) #that is, after the previous retrieval was completed
    if(debug){
      print(paste("retrieval number ",r))
    }


    current_retr <- retrievals[r,]
    scheduled_moment  <- chunks[chunks$name %in% current_retr$moment,]$created


    retrieval_moment$now <- pmax(available_moment$available_at, scheduled_moment )

    retrieval_cues <- current_retr[cue_names]

    history_till_now <- left_join(history,retrieval_moment,by= c("exp", "subj", "item"))
    #I remove the future words from the history
    history_till_now <- filter(history_till_now,accessed < now & !is.na(accessed) )
    
    ##########33
    ##retrieval using act-r formulas here:
    ret <- retrieve(cue_names,retrieval_cues, retrieval_moment,history_till_now)
    ##############

    if(any(is.na(ret$winner$accessed))){
      print(paste("retrieval",r))
      stop("accessed is = NA!!! at")
    }

    #I add the retrieval times to the history of events:
    history <- bind_rows(history,ret$winners)
    #arrange(history,exp,subj,item,accessed)

    #I remove the old now
    history<- select(history,-now)

    ret$actr_sim$retr <- current_retr$retr
    ret$actr_sim$retrieval_at <- current_retr$moment
    #arrange(ret$actr_sim,exp,subj,item)

 

    #this is when the retrieval is over:
    available_moment$available_at <- ret$winner$accessed

   
    #ASK FELIX IF I SHOULD LEAVE TIME FOR THE RULE
    
    result_list[[r]] <- ret$actr_sim
  }
    
    
  results<-bind_rows(result_list)   
 


 activations <- spread(
      select(results,exp,subj,item,chunkn,retr,retrieval_at,final_activation), 
    chunkn,final_activation)

 retrievals <- spread(
      select(results,exp,subj,item,chunkn,retr,retrieval_at,winner), 
    chunkn,winner)

  latencies <- summarise(
    group_by_(subset(results,winner),.dots=c("exp","subj","item","retrieval_at",actr_par_names))
    ,latency=sum(final_latency))


  complete_results =list(results=results,latencies=latencies,retrievals=retrievals,activations=activations, history=history)
  return(complete_results)
}



sim_actr <- function(conditions_list,actr_par=NULL,nsim=1,debug=F){
  if (is.null(actr_par)){
    print("Using default act-r parameters. See actr_default")
    actr_par <- actr_default_1subj
  }

  results <- plyr::llply(conditions_list, run_actr_each_condition, actr_par=actr_par,nsim=nsim,debug=debug)

  actr_sim <- list()
  actr_sim$results <- results
  actr_sim$actr_par <- actr_par
  actr_sim$nsim <- nsim

  #I don't know what's better
  #class(actr_sim) <- append(class(actr_sim),"actr")
  attr(actr_sim, "class") <- "actr"

  return(actr_sim)  
}



repeat_row <- function(data,...){
  data[ rep(seq_len(nrow(data)), ...),]  
}




actr_default_1subj <- as_data_frame(actr_defaults)





#' Print basic information regarding the actr simulation and a summary for the parameters of interest
#'
#' @param actr_sim an S3 actr object

summary.actr <- function(actr_sim){

  latencies <-  plyr::ldply(actr_sim$results, function(l){
    
    summarise(group_by(l$latencies,retrieval_at),Latency = mean(latency), SE= sd(latency)/sqrt(length(latency)))

  }) #end of plyr::llply


  summary <-  plyr::ldply(actr_sim$results, function(l){

     summarise(
      group_by(l$results,retr,chunkn,name,retrieval_at),Retr_P = mean(winner),
      Latency=mean(final_latency,na.rm=T),
      mean_activation=mean(final_activation))
  })
  
  retrievals <-  plyr::ldply(actr_sim$results, function(l){
      summarise_each(
      group_by(l$retrievals,retr,retrieval_at),funs(mean(.,na.rm=T)),matches("chunk|failure"))
  })
  
# spread(
#       select(results,exp,subj,item,chunkn,retr,retrieval_at,winner), 
#     chunkn,winner)



  nitem<-  plyr::llply(actr_sim$results, function(l){
    return(max(l$results$item))
  })


  lines <- c("# Summary of simulations",
    paste("# Number of simulations: ",actr_sim$nsim),
    paste(paste("# Number of items ",names(nitem)),nitem),
    paste("# Number of simulated subjects: ",nrow(actr_sim$actr_par),"\n"))

  cat(paste(lines,collapse="\n"))
  
  summary<- list(summary=summary, latencies=latencies,retrievals=retrievals)
  class(summary) <- append(class(summary),"sum_actr_sim")
  return(summary)
}


by_subj <- function(actr_sim,pars=NULL){

  latencies <-  plyr::ldply(actr_sim$results, function(l){
    
    summarise(group_by_(l$latencies,.dots=c("subj","retrieval_at",pars)),Latency = mean(latency), SE= sd(latency)/sqrt(length(latency)))

  }) 

  return(latencies)
}

abbreviate_UTF8 <- function(names.arg,minlength){
  abbreviate(iconv(names.arg, to='ASCII//TRANSLIT'),minlength=minlength,"both")
}

print.sum_actr_sim <- function(summary_actr, removeNaN=F, minlength = 10, ...) {
    if (minlength!=0) {
        summary_actr$summary <- mutate_each(summary_actr$summary,funs(abbreviate_UTF8(.,minlength)),-one_of("mean_activation","Retr_P","Latency"))
         
    }
      if(removeNaN){
           summary_actr$summary <- summary_actr$summary[!is.na(summary_actr$summary$Latency) & !is.nan(summary_actr$summary$Latency),]}

    print(unclass(summary_actr))
}




plot.actr <- function(actr_sim,pars=NULL){

  # by_subj <- ifelse(is.null(pars),FALSE,TRUE)
  # summary <- summary(actr_sim,latencies=TRUE,pars=pars[1],by_subj=by_subj) #only the first par

  if(is.null(pars)){
    summary <- summary(actr_sim)$latencies
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

    summary <- by_subj(actr_sim,pars=pars[1])
    if(length(pars)>1){warning("Using only the first parameter")}
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

 #more informative errors when something is wrong
sanity_check <- function(chunks,retrievals){


  #check for the obligatory columns in creation schedule
  must_cols_cs <- c("name","created","cat") 
  if(!all(must_cols_cs %in% colnames(chunks))){
     missing_cols <-must_cols_cs[which(!must_cols_cs %in% colnames(chunks))]
    stop(paste("Column",missing_cols,"missing from creation_schedule"))
  }

  #check for the obligatory columns in retrieval schedule
  must_cols_rs <- c("moment","cat") 
  if(!all(must_cols_rs %in% colnames(retrievals))){
     missing_cols <-must_cols_rs[which(!must_cols_rs %in% colnames(retrievals))]
    stop(paste("Column",missing_cols,"missing from retrievals_schedule"))
  }

  #check if the moments from the retrieval schedule exist in the creation schedule
  if(!all(retrievals$moment %in% chunks$name)){
    missing_name <-retrievals[which(!retrievals$moment %in% chunks$name),]$moment

     stop(paste("Moment '",missing_name,"' missing from creation_schedule's name"))
  }
 
}

