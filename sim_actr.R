library(plyr)

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

  schedule <- t(schedule)
  rownames(schedule) <- NULL
  schedule<-data.frame(schedule)

  
  #validation of the row names
  colnames(schedule) <- make.names(colnames(schedule))
  return(schedule)
  
}


run_actr_each_condition <- function(data_list,actr_par,nsim) {
  print("Starting to run model....")


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
  sim_trials <-  data.frame(
    item=rep(1:data_list$num_experimental_items,each=nrow(retrievals),times=nsim*nsubj),
    exp= rep(1:nsim,each=data_list$num_experimental_items*nrow(retrievals),times=nsubj ))
  #each retrieval for each subj for each item for each experiment:
  retr_subj <- actr_par[rep(seq_len(nrow(actr_par)), each=nrow(retrievals)*data_list$num_experimental_items*nsim),]
  subjects <- rep(seq_len(nrow(actr_par)), each=nrow(retrievals)*data_list$num_experimental_items*nsim)
  sim_trials <- cbind(sim_trials,retr_subj )
  sim_trials <- cbind(subj=subjects, sim_trials)
  rep_retr <- retrievals[rep(seq_len(nrow(retrievals)), data_list$num_experimental_items*nsim),]

  row.names(rep_retr) <- NULL
  row.names(sim_trials) <- NULL
  sim_retrievals <- cbind(rep_retr,sim_trials)


  #save here the histories of retrievals
  histories <- list()


  #this is the results of all the retrievals in all the simulations
  results <- ddply(sim_retrievals, .(subj,item,exp), function(a){
    #test with a <- sim_retrievals[sim_retrievals$item==1 & sim_retrievals$exp==1 & sim_retrievals$subj ==1,]
    #the first history for each itemxexp combination is just the creation schedule
    history <- data.frame(chunk=rownames(chunks),name=as.character(chunks$name),moment=chunks$created, stringsAsFactors = FALSE)  #I will add failure, so I need the columns to be characters and not factors

    #history:
  #   chunk name moment
  # 1 word1  NP1      0
  # 2 word2  DP2    400
  # 3 word3  VP2    600
  # 4 word4  NP2    800
  # 5 word5  NP3   1400
  # 6 word6  VP1   1800

    #moment of the first retrieval: look for retrieval moment when the item was created:
    
    available_moment <- -1 #before any retrieval
    res <- data.frame()
    for(r in 1:length(unique(a$retr))){
      #this checks if the retrieval can be done as schedule (a[r,]$moment) or delayed to the next available moment (available_moment) #that is, after the previous retrieval was completed
      scheduled_moment  <- chunks[chunks$name %in% a[r,]$moment,]$created

      retrieval_moment <- max(available_moment, scheduled_moment)
      
      actr_par_retr <- a[1,colnames(a) %in% colnames(actr_par)  ]
      retrieval_cues <- a[r,colnames(a) %in% cue_names]
      ret <- retrieve(cue_names,retrieval_cues , retrieval_moment,chunks,history,actr_par_retr )
      history <- ret$newhistory
      ret$output$retr <- as.character(a[r,]$retr)
      ret$output$retrieval_at <- as.character(a[r,]$moment)
      available_moment <- ret$finished_at

      #ASK FELIX IF I SHOULD LEAVE TIME FOR THE RULE
      res <- rbind(res, ret$output)
    }
    histories <<- c(histories,list(history))

    return(res)
  }) 

latencies <- ddply(subset(results,winner==1),.(exp,subj,item,retrieval_at),summarize,latency=sum(final_latency))

  complete_results =list(results=results,latencies=latencies,histories=histories)
  return(complete_results)
}


library(plyr)

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




summary.actr <- function(output,latencies=FALSE,by_subj=FALSE,removeNaN=FALSE){
    

 summ <-  ldply(output$results, function(l){
    df <- data.frame()
    if(latencies & by_subj){
    
      df<- ddply(l$latencies,.(subj,retrieval_at),summarize,mean_latency = mean(latency), SE= sd(latency)/sqrt(length(latency)))
    
    } else if(latencies & !by_subj){
      if(max(l$latencies$subj>1)){warning("SE is not corrected for subjects")}
      warning
      df <- ddply(l$latencies,.(retrieval_at),summarize,mean_latency = mean(latency), SE= sd(latency)/sqrt(length(latency)))
    
    } else if(!latencies & by_subj){
    
      df <- ddply(l$results,.(subj,retr,name,chunk),summarize,P=mean(winner),mean_latency=mean(final_latency,na.rm=T),mean_activation=mean(final_activation))
    
    } else if(!latencies & !by_subj){ 

      df<-  ddply(l$results,.(retr,name,chunk),summarize,P=mean(winner),mean_latency=mean(final_latency,na.rm=T),mean_activation=mean(final_activation))
      }
    return(df)
  })


  if(removeNaN){
    summ <- summ[!is.nan(summ$mean_latency),]
  }

  lines <- c("# Summary of simulations",
    paste("# Number of simulations: ",output$nsim),
    paste("# Number of simulated subjects: ",nrow(output$actr_par),"\n"))
  cat(paste(lines,collapse="\n"))
  return(summ)
}


plot.actr <- function(output){
  summary <- summary(output,latencies=TRUE)


  plot <- ggplot(summary, aes(x=.id, y=mean_latency, fill=.id))+ 
          facet_grid(. ~ retrieval_at) + 
          geom_bar(position=position_dodge(), stat="identity") +
          geom_errorbar(aes(ymin=mean_latency-2*SE, ymax=mean_latency+2*SE),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
          xlab("Experimental Condition") +
    ylab("Latency") +
    scale_fill_hue(name="Experimental Condition") +
    ggtitle("Summary") +
    #scale_y_continuous(breaks=seq(0:summary$)) +
    theme_bw()




  return(plot)
}