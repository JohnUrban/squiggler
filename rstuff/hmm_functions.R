## 
generate.statepath <- function(transitions, initial, states, length=10){
  statenums <- 1:length(states)
  current <- sample(statenums, size = 1, prob = initial)
  statepathseq <- c(states[current])
  length = length-1
  while(length > 0){
    upcoming <- sample(statenums, size=1, prob = transitions[current, ])
    current <- upcoming
    statepathseq <- c(statepathseq, states[current])
    length = length-1
  }
  return(paste0(statepathseq, collapse=""))
}

generate.emissions<- function(emissions, transitions, initial, states, length=10, withGaps=FALSE){
  if(!(withGaps)){## NOT gaps
    generate.emissions.fxn(emissions, transitions, initial, states, getseqfxn=get.sequence, length=length)
  } else {## gaps
    generate.emissions.fxn(emissions, transitions, initial, states, getseqfxn=get.sequence.withgaps, length=length)
  }
}

generate.emissions.fxn <- function(emissions, transitions, initial, states, getseqfxn=get.sequence, length=10){
  statenums <- 1:length(states)
  current <- sample(statenums, size = 1, prob = initial)
  statepathseq <- c(states[current])
  statepath <- c(current)
  emitted.data <- rnorm(n = 1, mean = emissions[1,current], sd = emissions[2,current])
  length = length-1
  while(length > 0){
    upcoming <- sample(statenums, size=1, prob = transitions[current, ])
    current <- upcoming
    statepathseq <- c(statepathseq, states[current])
    statepath <- c(statepath, current)
    emitted.data <- c(emitted.data, rnorm(n = 1, mean = emissions[1,current], sd = emissions[2,current]))
    length = length-1
  }#paste0(statepathseq, collapse="")
  seqinfo <- getseqfxn(states, statepath)
  return(list(statepath=statepath, statepathseq=seqinfo$seq, emitted.data=emitted.data, moves=seqinfo$moves))
}


generate.emissions.twoEmits <- function(emissions1, emissions2, transitions, initial, states, getseqfxn=get.sequence, length=10){
  statenums <- 1:length(states)
  current <- sample(statenums, size = 1, prob = initial)
  statepathseq <- c(states[current])
  statepath <- c(current)
  emitted.data1 <- rnorm(n = 1, mean = emissions1[1,current], sd = emissions1[2,current])
  emitted.data2 <- rnorm(n = 1, mean = emissions2[1,current], sd = emissions2[2,current])
  length = length-1
  while(length > 0){
    upcoming <- sample(statenums, size=1, prob = transitions[current, ])
    current <- upcoming
    statepathseq <- c(statepathseq, states[current])
    statepath <- c(statepath, current)
    emitted.data1 <- c(emitted.data1, rnorm(n = 1, mean = emissions1[1,current], sd = emissions1[2,current]))
    emitted.data2 <- c(emitted.data2, rnorm(n = 1, mean = emissions2[1,current], sd = emissions2[2,current]))
    length = length-1
  }#paste0(statepathseq, collapse="")
  seqinfo <- getseqfxn(states, statepath)
  return(list(statepath=statepath, statepathseq=seqinfo$seq, emitted.data1=emitted.data1, emitted.data2=emitted.data2, moves=seqinfo$moves))
}



forward <- function(emissions, transitions, initial, states, emitted.data){
  num.states <- length(states)
  num.emits <- length(emitted.data)
  Forward = matrix(data = rep(0, num.states*num.emits), nrow = num.states)
  
  ## initial
  for (i in 1:num.states){
    Forward[i, 1] <- initial[i]*dnorm(emitted.data[1], mean = emissions[1, i], sd = emissions[2, i])
  }
  
  ## iterate
  for(k in 2:num.emits){
    for (i in 1:num.states){
      for (m in 1:num.states){
        ##TODO
        #selection <- Forward[,j-1]*transitions[,i]
        fwd <- Forward[m,k-1]
        emit <- dnorm(emitted.data[k], mean = emissions[1, i], sd = emissions[2, i])
        trans <- transitions[m,i]
        Forward[i, k] <- Forward[i, k] + fwd*trans*emit
      }
    }
  }
  return(Forward)
}

# forward.long <- function(emissions, transitions, initial, states, emitted.data){
#   num.states <- length(states)
#   num.emits <- length(emitted.data)
#   Forward <- matrix(data = rep(0, num.states*num.emits), nrow = num.states)
#   scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
#   
#   ## initial
#   for (i in 1:num.states){
#     Forward[i, 1] <- initial[i]*dnorm(emitted.data[1], mean = emissions[1, i], sd = emissions[2, i])
#   }
#   ## scale to prevent underflow -- keep track of scaling
#   scalefactors[1,1] <- sum(Forward[, 1])
# #   scalefactors[2,1] <- scalefactors[1,1]
#   scalefactors[2,1] <- log(scalefactors[1,1])
#   Forward[,1] <- Forward[,1]/scalefactors[1,1]
# 
# #   ## iterate
# #   for(k in 2:num.emits){
# #     for (i in 1:num.states){
# #       for (m in 1:num.states){
# #         ##TODO
# #         #selection <- Forward[,j-1]*transitions[,i]
# #         fwd <- Forward[m,k-1]
# #         emit <- dnorm(emitted.data[k], mean = emissions[1, i], sd = emissions[2, i])
# #         trans <- transitions[m,i]
# #         Forward[i, k] <- Forward[i, k] + fwd*trans*emit
# #       }
# #     }
# #     scalefactors[1,k] <- sum(Forward[, k])
# # #     scalefactors[2,k] <- scalefactors[1,k]*scalefactors[2,k-1]
# #     scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k-1]
# #     Forward[,k] <- Forward[,k]/scalefactors[1,k]
# #   }
# 
#    ## iterate -- partially vectorized (m loop gone, I have a new version where both i and m loops are gone)
#   for(k in 2:num.emits){
#     for (i in 1:num.states){
#       emit <- dnorm(emitted.data[k], mean = emissions[1, i], sd = emissions[2, i])
#       trans <- trans <- transitions[,i]
#       fwd <- Forward[,k-1]
#       Forward[i, k] <- emit*sum(fwd*trans)
#     }
#     scalefactors[1,k] <- sum(Forward[, k])
#     scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k-1]
#     Forward[,k] <- Forward[,k]/scalefactors[1,k]
#   }
# 
#   return(list(forward=Forward, scales=scalefactors))
#   ## mutiply forward column by row2,samecol in scale factors OR by exp(row3,samecol in scalfactors) to get actual value for forward
#   ## 
# }


# forward.long.twoEmits <- function(emissions1, emissions2, transitions, initial, states, emitted.data1, emitted.data2){
#   num.states <- length(states)
#   num.emits <- length(emitted.data1)
#   Forward <- matrix(data = rep(0, num.states*num.emits), nrow = num.states)
#   scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
#   
#   ## initial
#   for (i in 1:num.states){
#     ep1 <- dnorm(emitted.data1[1], mean = emissions1[1, i], sd = emissions1[2, i])
#     ep2 <- dnorm(emitted.data2[1], mean = emissions2[1, i], sd = emissions2[2, i])
#     Forward[i, 1] <- initial[i]*ep1*ep2
#   }
#   ## scale to prevent underflow -- keep track of scaling
#   scalefactors[1,1] <- sum(Forward[, 1])
#   #   scalefactors[2,1] <- scalefactors[1,1]
#   scalefactors[2,1] <- log(scalefactors[1,1])
#   Forward[,1] <- Forward[,1]/scalefactors[1,1]
#   
#   ## iterate
#   for(k in 2:num.emits){
#     for (i in 1:num.states){
#       for (m in 1:num.states){
#         ##TODO
#         #selection <- Forward[,j-1]*transitions[,i]
#         fwd <- Forward[m,k-1]
#         ep1 <- dnorm(emitted.data1[k], mean = emissions1[1, i], sd = emissions1[2, i])
#         ep2 <- dnorm(emitted.data2[k], mean = emissions2[1, i], sd = emissions2[2, i])
#         trans <- transitions[m,i]
#         Forward[i, k] <- Forward[i, k] + fwd*trans*ep1*ep2
#       }
#     }
#     scalefactors[1,k] <- sum(Forward[, k])
#     #     scalefactors[2,k] <- scalefactors[1,k]*scalefactors[2,k-1]
#     scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k-1]
#     Forward[,k] <- Forward[,k]/scalefactors[1,k]
#   }
#   return(list(forward=Forward, scales=scalefactors))
#   ## mutiply forward column by row2,samecol in scale factors OR by exp(row3,samecol in scalfactors) to get actual value for forward
#   ## 
# }

backward <- function(emissions, transitions, initial, states, emitted.data){
  num.states <- length(states)
  num.emits <- length(emitted.data)
  Backward = matrix(data = rep(0, num.states*num.emits), nrow = num.states)
  
  ## initial
  for (i in 1:num.states){
    Backward[i, num.emits] <- 1
  }
  
  ## iterate
  for(k in (num.emits-1):1){
    for (i in 1:num.states){
      del<-0
      for (m in 1:num.states){
        back = Backward[m, k+1]
        emit = dnorm(emitted.data[k+1], mean = emissions[1, m], sd = emissions[2, m])
        trans = transitions[i,m]
        Backward[i,k] <- Backward[i,k] + trans*emit*back
      }  
    }
  }
  return(Backward)
}


# backward.long <- function(emissions, transitions, initial, states, emitted.data){
#   num.states <- length(states)
#   num.emits <- length(emitted.data)
#   Backward = matrix(data = rep(0, num.states*num.emits), nrow = num.states)
#   scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
#   
#   ## initial
#   for (i in 1:num.states){
#     Backward[i, num.emits] <- 1
#   }
#   ## scale to prevent underflow -- keep track of scaling
#   scalefactors[1,num.emits] <- sum(Backward[, num.emits])
# #   scalefactors[2,num.emits] <- scalefactors[1,num.emits]
#   scalefactors[2,num.emits] <- log(scalefactors[1,num.emits])
#   Backward[,num.emits] <- Backward[,num.emits]/scalefactors[1,num.emits]
#   
#   ## iterate
#   for(k in (num.emits-1):1){
#     for (i in 1:num.states){
#       del<-0
#       for (m in 1:num.states){
#         back = Backward[m, k+1]
#         emit = dnorm(emitted.data[k+1], mean = emissions[1, m], sd = emissions[2, m])
#         trans = transitions[i,m]
#         Backward[i,k] <- Backward[i,k] + trans*emit*back
#       }  
#     }
#     scalefactors[1,k] <- sum(Backward[, k])
# #     scalefactors[2,k] <- scalefactors[1,k]*scalefactors[2,k+1]
#     scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k+1]
#     Backward[,k] <- Backward[,k]/scalefactors[1,k]
#   }
#   return(list(backward=Backward, scales=scalefactors))
# }
##VECTORIZED BACKWARD.LONG BELOW


# backward.long.twoEmits <- function(emissions1, emissions2, transitions, initial, states, emitted.data1, emitted.data2){
#   num.states <- length(states)
#   num.emits <- length(emitted.data1)
#   Backward = matrix(data = rep(0, num.states*num.emits), nrow = num.states)
#   scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
#   
#   ## initial
#   for (i in 1:num.states){
#     Backward[i, num.emits] <- 1
#   }
#   ## scale to prevent underflow -- keep track of scaling
#   scalefactors[1,num.emits] <- sum(Backward[, num.emits])
#   #   scalefactors[2,num.emits] <- scalefactors[1,num.emits]
#   scalefactors[2,num.emits] <- log(scalefactors[1,num.emits])
#   Backward[,num.emits] <- Backward[,num.emits]/scalefactors[1,num.emits]
#   
#   ## iterate
#   for(k in (num.emits-1):1){
#     for (i in 1:num.states){
#       del<-0
#       for (m in 1:num.states){
#         back = Backward[m, k+1]
#         ep1 = dnorm(emitted.data1[k+1], mean = emissions1[1, m], sd = emissions1[2, m])
#         ep2 = dnorm(emitted.data2[k+1], mean = emissions2[1, m], sd = emissions2[2, m])
#         trans = transitions[i,m]
#         Backward[i,k] <- Backward[i,k] + trans*ep1*ep2*back
#       }  
#     }
#     scalefactors[1,k] <- sum(Backward[, k])
#     #     scalefactors[2,k] <- scalefactors[1,k]*scalefactors[2,k+1]
#     scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k+1]
#     Backward[,k] <- Backward[,k]/scalefactors[1,k]
#   }
#   return(list(backward=Backward, scales=scalefactors))
# }

####### UPDATED/COMPLETELY.VECTORIZED FORWARD AND BACKWARD ALGOS ####################
## new fwd.long pulls emit calculation out side of m loop b/c it only depends on i and k and needs tobe calc just once
## This needs to be multiplied by every t_m->i*f_m,k-1
forward.long <- function(emissions, transitions, initial, states, emitted.data){
  num.states <- length(states)
  num.emits <- length(emitted.data)
  Forward <- matrix(data = rep(0, num.states*num.emits), nrow = num.states)
  scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
  
  ## initial
  Forward[, 1] <- initial*dnorm(emitted.data[1], mean = emissions[1, ], sd = emissions[2, ])
  
  ## scale to prevent underflow -- keep track of scaling
  scalefactors[1,1] <- sum(Forward[, 1])
  scalefactors[2,1] <- log(scalefactors[1,1])
  Forward[,1] <- Forward[,1]/scalefactors[1,1]

  ## iterate
  for(k in 2:num.emits){
    emit <- dnorm(emitted.data[k], mean = emissions[1, ], sd = emissions[2, ])
    Forward[, k] <- emit* Forward[,k-1] %*% transitions ## same as emit* Forward[,k-1] * colSums(transitions)
#     a <- Forward[,k-1] %*% transitions
# Forward[, k] <- emit*a
# #     print(emit)
# #     print(emit*Forward[,k-1])
# #     print(Forward[,k])
# #     print("")
    scalefactors[1,k] <- sum(Forward[, k])
    scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k-1]
    Forward[,k] <- Forward[,k]/scalefactors[1,k]
  }
  
  return(list(forward=Forward, scales=scalefactors))
  ## mutiply forward column by row2,samecol in scale factors OR by exp(row3,samecol in scalfactors) to get actual value for forward
  ## 
}

forward.long.twoEmits <- function(emissions1, emissions2, transitions, initial, states, emitted.data1, emitted.data2){
  num.states <- length(states)
  num.emits <- length(emitted.data1)
  Forward <- matrix(data = rep(0, num.states*num.emits), nrow = num.states)
  scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
  
  ## initial
  ep1 <- dnorm(emitted.data1[1], mean = emissions1[1, ], sd = emissions1[2, ])
  ep2 <- dnorm(emitted.data2[1], mean = emissions2[1, ], sd = emissions2[2, ])
  Forward[, 1] <- initial*ep1*ep2
  
  ## scale to prevent underflow -- keep track of scaling
  scalefactors[1,1] <- sum(Forward[, 1])
  scalefactors[2,1] <- log(scalefactors[1,1])
  Forward[,1] <- Forward[,1]/scalefactors[1,1]
  
  ## iterate
  for(k in 2:num.emits){
    emit <- dnorm(emitted.data[k], mean = emissions[1, ], sd = emissions[2, ])
    ep1 <- dnorm(emitted.data1[k], mean = emissions1[1, ], sd = emissions1[2, ])
    ep2 <- dnorm(emitted.data2[k], mean = emissions2[1, ], sd = emissions2[2, ])
    Forward[, k] <- ep1*ep2* Forward[,k-1] %*% transitions ## same as emit* Forward[,k-1] * colSums(transitions)
    scalefactors[1,k] <- sum(Forward[, k])
    scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k-1]
    Forward[,k] <- Forward[,k]/scalefactors[1,k]   
    scalefactors[1,k] <- sum(Forward[, k])
    scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k-1]
    Forward[,k] <- Forward[,k]/scalefactors[1,k]
  }
  return(list(forward=Forward, scales=scalefactors))
  ## mutiply forward column by row2,samecol in scale factors OR by exp(row3,samecol in scalfactors) to get actual value for forward
  ## 
}

#### VECTORIZED BACKWARDS ########
backward.long <- function(emissions, transitions, initial, states, emitted.data){
  num.states <- length(states)
  num.emits <- length(emitted.data)
  Backward = matrix(data = rep(0, num.states*num.emits), nrow = num.states)
  scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
  
  ## initial
  Backward[ , num.emits] <- 1
  
  ## scale to prevent underflow -- keep track of scaling
  scalefactors[1,num.emits] <- sum(Backward[, num.emits])
  scalefactors[2,num.emits] <- log(scalefactors[1,num.emits])
  Backward[,num.emits] <- Backward[,num.emits]/scalefactors[1,num.emits]
  
  ## iterate
  for(k in (num.emits-1):1){
    emit <- matrix(dnorm(emitted.data[k+1], mean = emissions[1, ], sd = emissions[2, ]))
#     print(Backward[, k+1] * emit)
    Backward [, k] <- transitions %*% (Backward[, k+1] * emit)
    scalefactors[1,k] <- sum(Backward[, k])
    scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k+1]
    Backward[,k] <- Backward[,k]/scalefactors[1,k]
  }
  return(list(backward=Backward, scales=scalefactors))
}

prob.data <- function(Forward){
  num.emits <- dim(Forward)[2]
  return(sum(Forward[,num.emits]))
}

backward.long.twoEmits <- function(emissions1, emissions2, transitions, initial, states, emitted.data1, emitted.data2){
  num.states <- length(states)
  num.emits <- length(emitted.data1)
  Backward = matrix(data = rep(0, num.states*num.emits), nrow = num.states)
  scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
  
  ## initial
  Backward[ , num.emits] <- 1
  
  ## scale to prevent underflow -- keep track of scaling
  scalefactors[1,num.emits] <- sum(Backward[, num.emits])
  scalefactors[2,num.emits] <- log(scalefactors[1,num.emits])
  Backward[,num.emits] <- Backward[,num.emits]/scalefactors[1,num.emits]
  
  ## iterate
  for(k in (num.emits-1):1){
    ep1 = dnorm(emitted.data1[k+1], mean = emissions1[1, ], sd = emissions1[2, ])
    ep2 = dnorm(emitted.data2[k+1], mean = emissions2[1, ], sd = emissions2[2, ])
    Backward [, k] <- transitions %*% (Backward[, k+1] * ep1 * ep2)
    scalefactors[1,k] <- sum(Backward[, k])
    scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k+1]
    Backward[,k] <- Backward[,k]/scalefactors[1,k]
    scalefactors[1,k] <- sum(Backward[, k])
    scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k+1]
    Backward[,k] <- Backward[,k]/scalefactors[1,k]
  }
  return(list(backward=Backward, scales=scalefactors))
}




############## PROB DATA ##########################

prob.data.long <- function(Forward, scalefactors){
  num.emits <- dim(Forward)[2]
#   return(sum(Forward[,num.emits])*scalefactors[2,num.emits]) 
  return(sum(Forward[,num.emits])*exp(scalefactors[2,num.emits]))
}


# viterbi <- function(emissions, transitions, initial, states, emitted.data){
#   num.states <- length(states)
#   num.emits <- length(emitted.data)
#   pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
#   Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)   
#   ## need to add log_probs instead of multiply probs to prevent underflow
# 
#   for (i in 1:num.states){
#     Viterbi[i,1] <- log(initial[i]) + dnorm(emitted.data[1], mean = emissions[1, i], sd = emissions[2, i], log = TRUE)
#     pointer[1,i] <- 1
#   }
# 
#   for (j in 2:num.emits){
#     for (i in 1:num.states){
#       selection <- Viterbi[,j-1] + log(transitions[,i])
#       maxstate <- which.max(selection)
#       Viterbi[i,j] <- dnorm(emitted.data[j], mean = emissions[1, i], sd = emissions[2, i], log = TRUE) + selection[maxstate]
#       pointer[j,i] <- maxstate 
#     }
#   }     
#   
#   viterbi_path <- matrix(data = rep(0, num.emits), nrow = 1)
#   viterbi_path[num.emits] <- which.max(Viterbi[,num.emits]);  
#   viterbi_prob <- Viterbi[viterbi_path[num.emits], num.emits]
# 
#   for (j in num.emits:2){
#     viterbi_path[j-1] <- pointer[j,viterbi_path[j]]
#   }
#   return(list(viterbi_path=viterbi_path, viterbi_seq=get.sequence(states, viterbi_path), viterbi_prob=viterbi_prob))
# }

## new viterbi code, turns all probs into logprobs just once before iterating instead of logging every time
# viterbi <- function(emissions, transitions, initial, states, emitted.data, getseqfxn=get.sequence, logprobs=FALSE){
#   ## logprobs = whether or not transition and initial matrices are logged yet, assumes false
#   if(!(logprobs)){
#     initial = log(initial)
#     transitions = log(transitions)
#   }
#   num.states <- length(states)
#   num.emits <- length(emitted.data)
#   pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
#   Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)   
#   ## need to add log_probs instead of multiply probs to prevent underflow
# 
#   for (i in 1:num.states){
#     Viterbi[i,1] <- initial[i] + dnorm(emitted.data[1], mean = emissions[1, i], sd = emissions[2, i], log = TRUE)
#     pointer[1,i] <- 1
#   }
#   
#   for (j in 2:num.emits){
#     for (i in 1:num.states){
#       selection <- Viterbi[,j-1] + transitions[,i]
#       maxstate <- which.max(selection)
#       Viterbi[i,j] <- dnorm(emitted.data[j], mean = emissions[1, i], sd = emissions[2, i], log = TRUE) + selection[maxstate]
#       pointer[j,i] <- maxstate 
#     }
#   }     
#   
#   viterbi_path <- matrix(data = rep(0, num.emits), nrow = 1)
#   viterbi_path[num.emits] <- which.max(Viterbi[,num.emits]);  
#   viterbi_prob <- Viterbi[viterbi_path[num.emits], num.emits]
#   
#   for (j in num.emits:2){
#     viterbi_path[j-1] <- pointer[j,viterbi_path[j]]
#   }
#   seqinfo=getseqfxn(states, viterbi_path)
#   return(list(viterbi_path=viterbi_path, viterbi_seq=seqinfo$seq, viterbi_prob=viterbi_prob, moves=seqinfo$moves))
# }

##VECTORIZED VITERBI -- actually slower....
viterbi <- function(emissions, transitions, initial, states, emitted.data, getseqfxn=get.sequence, logprobs=FALSE){
  ## logprobs = whether or not transition and initial matrices are logged yet, assumes false
  if(!(logprobs)){
    initial = log(initial)
    transitions = log(transitions)
  }
  num.states <- length(states)
  num.emits <- length(emitted.data)
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)   
  ## need to add log_probs instead of multiply probs to prevent underflow
  
  Viterbi[ ,1] <- initial + dnorm(emitted.data[1], mean = emissions[1, ], sd = emissions[2, ], log = TRUE)
  pointer[1, ] <- 1
  f <- function(x){i <- which.max(x); y <- x[i]; return(c(i,y))}
  for (j in 2:num.emits){
    selection <- Viterbi[,j-1] + transitions
#     maxstates <- apply(selection, 2, f)
#     Viterbi[ ,j] <- dnorm(emitted.data[j], mean = emissions[1, ], sd = emissions[2, ], log = TRUE) + maxstates[2,]
#     pointer[j, ] <- maxstates[1,]
    for (i in 1:num.states){
#       selection <- Viterbi[,j-1] + transitions[,i]
#       maxstate <- which.max(selection)
      maxstate <- which.max(selection[,i])
      Viterbi[i,j] <- dnorm(emitted.data[j], mean = emissions[1, i], sd = emissions[2, i], log = TRUE) + selection[maxstate,i]
      pointer[j,i] <- maxstate 
    }
  }     
  
  viterbi_path <- matrix(data = rep(0, num.emits), nrow = 1)
  viterbi_path[num.emits] <- which.max(Viterbi[,num.emits]);  
  viterbi_prob <- Viterbi[viterbi_path[num.emits], num.emits]
  
  for (j in num.emits:2){
    viterbi_path[j-1] <- pointer[j,viterbi_path[j]]
  }
  seqinfo=getseqfxn(states, viterbi_path)
  return(list(viterbi_path=viterbi_path, viterbi_seq=seqinfo$seq, viterbi_prob=viterbi_prob, moves=seqinfo$moves))
}


viterbi.twoEmits <- function(emissions1, emissions2, transitions, initial, states, emitted.data1, emitted.data2, getseqfxn=get.sequence){
  num.states <- length(states)
  num.emits <- length(emitted.data1)
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)   
  ## need to add log_probs instead of multiply probs to prevent underflow
  
  for (i in 1:num.states){
    ep1 <- dnorm(emitted.data1[1], mean = emissions1[1, i], sd = emissions1[2, i], log = TRUE)
    ep2 <- dnorm(emitted.data2[1], mean = emissions2[1, i], sd = emissions2[2, i], log = TRUE)
    Viterbi[i,1] <- log(initial[i]) + ep1 + ep2
    pointer[1,i] <- 1
  }
  
  for (j in 2:num.emits){
    for (i in 1:num.states){
      selection <- Viterbi[,j-1] + log(transitions[,i])
      maxstate <- which.max(selection)
      ep1 <- dnorm(emitted.data1[j], mean = emissions1[1, i], sd = emissions1[2, i], log = TRUE)
      ep2 <- dnorm(emitted.data2[j], mean = emissions2[1, i], sd = emissions2[2, i], log = TRUE)
      Viterbi[i,j] <- ep1 + ep2 + selection[maxstate]
      pointer[j,i] <- maxstate 
    }
  }     
  
  viterbi_path <- matrix(data = rep(0, num.emits), nrow = 1)
  viterbi_path[num.emits] <- which.max(Viterbi[,num.emits]);  
  viterbi_prob <- Viterbi[viterbi_path[num.emits], num.emits]
  
  for (j in num.emits:2){
    viterbi_path[j-1] <- pointer[j,viterbi_path[j]]
  }
  seqinfo=getseqfxn(states, viterbi_path)
  return(list(viterbi_path=viterbi_path, viterbi_seq=seqinfo$seq, viterbi_prob=viterbi_prob, moves=seqinfo$moves))
}


compare.statepath <- function(sp1, sp2){
  # where sp is a numeric vector or char vector with each state its own element -- not seq format
  total <- length(sp1)
  ident <- sum(sp1 == sp2)
  edit.dist <- total - ident
  return(list(edit.dist=edit.dist, identical.count=ident, pct.id=100*ident/total))
}

compare.seq <- function(s1, s2){
  ## assumes seqs same len
  un.ident=0
  un.total <- nchar(s1)
  for(i in 1:un.total){un.ident <- (un.ident + (substr(s1, i, i) == substr(s2,i,i)))}
  unaligned.edit.dist <- un.total - un.ident
  # align
  sigma <- nucleotideSubstitutionMatrix(match = 0, mismatch = -1, baseOnly = TRUE)
  x <- pairwiseAlignment(s1, s2, substitutionMatrix = sigma, gapOpening = -1, gapExtension = -1, scoreOnly = FALSE)
  ident <- 0
  total <- nchar(subject(x))
  for(i in 1:total){ident <- (ident + (substr(subject(x), i, i) == substr(pattern(x),i,i)))}
  edit.dist <- total - ident
  return(list(unaligned.edit.dist=unaligned.edit.dist, unaligned.identical.count=un.ident, unaligned.pct.id=100*un.ident/un.total, aln.edit.dist=edit.dist, aln.identical.count=ident, aln.pct.total=100*ident/total, aln=x))
}

centroid <- function(Forward, Backward, states, getseqfxn=get.sequence){
  ## F and B matrices are from small sequence fxns -- not scaled
  num.states <- length(states)
  num.emits <- dim(Forward)[2]
  centroid.path <- matrix(data = rep(0, num.emits), nrow = 1)
  probs <- matrix(data = rep(0, num.emits), nrow = 1)
  pd <- prob.data(Forward = Forward)
  for (i in 1:num.emits){
    fb <- Forward[,i]*Backward[,i]
    max.state <- which.max(fb)
    centroid.path[i] <- max.state
#     probs[i] <- max(fb)/pd ## should be divided by prob.data...?
  }
  seqinfo = getseqfxn(states, centroid.path)
  return(list(centroid.path=centroid.path, probs=probs, centroid.seq=seqinfo$seq, moves=seqinfo$moves))

}


centroid.long <- function(Forward, F_scales, Backward, B_scales, states, getseqfxn=get.sequence){
  ##F and B are scaled long seq matrices -- the scales are scalefactors that come with them out of long fxns
  num.states <- length(states)
  num.emits <- dim(Forward)[2]
  centroid.path <- matrix(data = rep(0, num.emits), nrow = 1)
  logprobs <- matrix(data = rep(0, num.emits), nrow = 1)
  string <- c()
  for (i in 1:num.emits){
    fb <- Forward[,i]*Backward[,i]
    max.state <- which.max(fb)
    centroid.path[i] <- max.state
    #     logprobs[i] <- exp(F_scales[i])*exp(B_scales[i])*max(fb) ## this needs scaling -- # should be divided by prob.data...?
    ## probs for short and long do not yet agree ##TODO
  }
  seqinfo = getseqfxn(states, centroid.path, checkoverlap=FALSE, posterior.decoding=TRUE)
  #   print(states)
  #   print(centroid.path)
  return(list(centroid.path=centroid.path, logprobs=logprobs, centroid.seq=seqinfo$seq, moves=seqinfo$moves))
  # posterior decoding allows most probable kmers in each position to be non-overlapping
  # however, it turns out it gets it mostly right, so if you just build a sequence from the most probable kmers
  #   it still only gets the bases added to the sequence wrong <= num states it got wrong
  #   however, it also turns out (w/ viterbi too) that the sequence is more accurate than the state path
  #   b/c when a state is wrong, the sequence still has a 25% chance of the correct base being there anyway
  # can also just use kmers and probs of each state to compare to viterbi -- statepath still possible to compare
  # and state path still does better overall interestingly... its just that right now I have not figured out a way
  # to turn the state path into a DNA sequence that is also better than the viterbi sequence.
  # the reason I had it working at all is because I was just pasting kmers together
  # 1 solution may be to just not check for overlap and build the sequence
  
}





simulateACGT <- function(n=100, seqlen=5000){
  cwins <- 0
  vwins <- 0
  for (i in 1:n){
    emissions = data.frame(A=c(rnorm(1,67,3),rnorm(1,1,0.1)), C=c(rnorm(1,67,3),rnorm(1,1,0.1)), G=c(rnorm(1,67,3),rnorm(1,1,0.1)), T=c(rnorm(1,67,3),rnorm(1,1,0.1)))
    A = rpois(4,10); A = A/sum(A)
    C = rpois(4,10); C = C/sum(C)
    G = rpois(4,10); G = G/sum(G)
    T = rpois(4,10); T = T/sum(T)
    transitions = matrix(data = c(A,C,G,T), nrow = 4, byrow = TRUE)
    initial = rpois(4,10); initial = initial/sum(initial)
    states = c("A","C","G","T")
    answer <- generate.emissions(emissions, transitions, initial, states, length=seqlen) 
    fl <- forward.long(emissions, transitions, initial, states, answer$emitted.data)
    dl <- prob.data.long(Forward = fl$forward, scalefactors = fl$scales)
    bl <- backward.long(emissions, transitions, initial, states, answer$emitted.data)
    v <- viterbi(emissions, transitions, initial, states, answer$emitted.data)
    cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales)
    vscore <- compare.statepath(v$viterbi_path, answer$statepath) 
    cscore <- compare.statepath(cl$centroid.path, answer$statepath)  
    if (vscore$pct.id > cscore$pct.id){vwins <- vwins+1}
    else {cwins <- cwins + 1}
  }
  return(list(centrod.wins=cwins, vertibi.wins=vwins))
}


get.sequence <- function(states, statepath, checkoverlap=TRUE, posterior.decoding=TRUE){
  ## states are some type of kmer
  ## statepath is vector of numbers (indexes)
  moves = rep(1, length(statepath)) ## first move is 0
  moves[1] <- 0
  k <- nchar(states[1])
  if(k == 1){return(list(seq=paste0(states[statepath],collapse=""), moves=moves))}
  else{
    #init
    seq <- c(states[statepath[1]])
    ##iter
    for(i in 2:length(statepath)){
      lastSuffix <- substr(states[statepath[i-1]],2,k)
      currentPrefix <- substr(states[statepath[i]],1,k-1)
      if(lastSuffix != currentPrefix && checkoverlap){return(c(paste(c("Error: suffix != prefix on ", i, "th move.")), seq))}
      seq <- c(seq, substr(states[statepath[i]],k,k))
#       moves[i] <- 1
    }
    return(list(seq=paste0(seq, collapse=""), moves=moves))
  }
}


generate.kmer.transition.probs <- function(states, nonzero.trans=NA){
  ## can set nonzero.trans to any vector -- for DNA length 4
  k <- nchar(states[1])
  if(k == 1){## assumes num.states is 4 -- may have to change
    if (!(is.na(nonzero.trans))){
      transitions = matrix(data = rep(nonzero.trans, 4), nrow = 4, byrow = TRUE) 
    }else{
      A = rpois(4,10); A = A/sum(A)
      C = rpois(4,10); C = C/sum(C)
      G = rpois(4,10); G = G/sum(G)
      T = rpois(4,10); T = T/sum(T)
      transitions = matrix(data = c(A,C,G,T), nrow = 4, byrow = TRUE) 
    }
    return(transitions)
  }
  else{
    num.states <- length(states)
    transitions <- matrix(rep(0, num.states*num.states), nrow = num.states)
    prefix <- list()
    for (i in 1:num.states){
      pref <- substr(states[i], 1,k-1)
      prefix[[pref]] <- c(prefix[[pref]],i)
    }

    if(!(is.na(nonzero.trans))){
      nonzero.trans <- nonzero.trans/sum(nonzero.trans)
    }
    for(i in 1:num.states){
      current.suffix <- substr(states[i], 2, k)
      if(is.na(nonzero.trans)){
        trans <- rpois(4,10)
        trans <- trans/sum(trans)
      }
      t <- 1 
      for (j in prefix[[current.suffix]]){
        if(is.na(nonzero.trans)){
          transitions[i,j] <- trans[t]
        }else{
          transitions[i,j] <- nonzero.trans[t]
        }
        t <- t+1 
      }
    }
    return(transitions)
  }
}

generate.kmer.emission.probs <- function(states, level=TRUE){
  # mu.mean and sigma.mean are the mean and std dev of the r7.3 state level means to be used to generate emission means 
  # mu.sd, sigma.sd -- same for std devs of signals
  
  if(level){
    mu.mean <- 65.57454
    sigma.mean <- 6.497453
    mu.sd <- 1.163836
    sigma.sd <- 0.4116285
  } else { ## sd emission
    mu.mean <- 1.37316
    sigma.mean <- 0.3144043
    mu.sd <- 0.1761904
    sigma.sd <- 0.06263217
  }
  ####
  num.states <- length(states)
  emissions <- matrix(rep(0, num.states*2), nrow=2)
  for (i in 1:num.states){
    emissions[1,i] <- rnorm(1,mu.mean,sigma.mean)
    emissions[2,i] <- abs(rnorm(1,mu.sd,sigma.sd))
  }
  return(emissions)
}

########################################################################
## TODO: transitions including staying in place and skipping

generate.kmer.transition.probs.withgaps <- function(states){
  ## can set nonzero.trans to any vector -- for DNA length 4
  k <- nchar(states[1])
  if(k == 1 || k == 2){## actually dimers should be allowed to visit any other dimer for gaps (e.g. if its not one of the 4 dimers that ovelap current dimer)
#     return(generate.kmer.transition.probs(states))
    return("This has been designed to work with three-mers as the smallest kmer. \n")
  }
  else{
    num.states <- length(states)
    transitions <- matrix(rep(0, num.states*num.states), nrow = num.states)
    prefix <- list()
    ## 1 to k-1 and k-2
    for (i in 1:num.states){
      pref <- substr(states[i], 1,k-1)
      prefix[[pref]] <- c(prefix[[pref]],i)
      pref <- substr(states[i], 1,k-2)
      prefix[[pref]] <- c(prefix[[pref]],i)
    }

    ## create
    for(i in 1:num.states){
      ## overlap by 4 (move = 1)
      current.suffix <- substr(states[i], 2, k)
      trans <- rpois(4,365)
      #       trans <- trans/sum(trans)
      t <- 1 
      for (j in prefix[[current.suffix]]){
        transitions[i,j] <- trans[t]
        t <- t+1 
      }
      ## overlap by 3 (move = 2) -- add additional counts
      current.suffix <- substr(states[i], 3, k)
      trans <- rpois(16,4)
      t <- 1
      for (j in prefix[[current.suffix]]){
        transitions[i,j] <- transitions[i,j] + trans[t]
        t <- t+1 
      }
      ## add additional probability to staying in place (move = 0)
      current.suffix <- states[i]
      trans <- rpois(1,20)
      transitions[i,i] <- transitions[i,i] + trans
      ## normalize all counts by sum to create probs that sum to 1
      transitions[i,] <- transitions[i,]/sum(transitions[i,])
    }
    return(transitions)
  }
}



get.sequence.withgaps <- function(states, statepath, checkoverlap=TRUE, posterior.decoding=FALSE){##TODO -- only copied get.sequence so far
  ## states are some type of kmer
  ## statepath is vector of numbers (indexes)
  moves = rep(0, length(statepath)) ## first move is 0
  k <- nchar(states[1])
  if(k == 1 || k == 2){
#     return(paste0(states[statepath],collapse=""))
    return("This currently only works with 3-mers as smallest kmer.")
  }
  else{
    #init
    seq <- c(states[statepath[1]])
    moves[1] <- 0
    #iter
    for(i in 2:length(statepath)){
      lastSuffix <- substr(states[statepath[i-1]],2,k)
      currentPrefix <- substr(states[statepath[i]],1,k-1)
      if(lastSuffix == currentPrefix){
        seq <- c(seq, substr(states[statepath[i]],k,k))
        moves[i] <- 1
      } else {
        lastSuffix <- substr(states[statepath[i-1]],3,k)
        currentPrefix <- substr(states[statepath[i]],1,k-2)
        if(lastSuffix == currentPrefix){
          seq <- c(seq, substr(states[statepath[i]],k-1,k))
          moves[i] <- 2
        } else if (statepath[i-1] == statepath[i]) { 
          ## by checking same state last, only heteropolymers affected
          ## homopolymers would be caught in first condition
          moves[i] <- 0
          ## nothing is added to sequence
          ## could make another fxn that just spits out events and states line by line like 'template events' in f5
        } 
        ## ELSE::: do what? ... in other one just added centroid seq regardless...
        else if (posterior.decoding){
          seq <- c(seq, substr(states[statepath[i]],k,k))
          moves[i] <- -1 
          ## -1 means it was an "illegal" move (move to a kmer that does not overlap by k-1 or k-2)
          ## it turns out that adding the base from the illegal move does not hurt the seq overall much
        }
      }
      
    }
    return(list(seq=paste0(seq, collapse=""), moves=moves))
  }
}



generate.kmer.initial.probs <- function(states, uniform=TRUE){
  num.states <- length(states)
  if(uniform){
    initial<- rep(1/num.states, num.states)
  }else{
    initial <- rpois(num.states, 10)
    initial <- initial/sum(initial)
  }
  return(initial)
}



dimers <-c('AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT')

trimers <- c('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT')

fourmers <- c('AAAA', 'AAAC', 'AAAG', 'AAAT', 'AACA', 'AACC', 'AACG', 'AACT', 'AAGA', 'AAGC', 'AAGG', 'AAGT', 'AATA', 'AATC', 'AATG', 'AATT', 'ACAA', 'ACAC', 'ACAG', 'ACAT', 'ACCA', 'ACCC', 'ACCG', 'ACCT', 'ACGA', 'ACGC', 'ACGG', 'ACGT', 'ACTA', 'ACTC', 'ACTG', 'ACTT', 'AGAA', 'AGAC', 'AGAG', 'AGAT', 'AGCA', 'AGCC', 'AGCG', 'AGCT', 'AGGA', 'AGGC', 'AGGG', 'AGGT', 'AGTA', 'AGTC', 'AGTG', 'AGTT', 'ATAA', 'ATAC', 'ATAG', 'ATAT', 'ATCA', 'ATCC', 'ATCG', 'ATCT', 'ATGA', 'ATGC', 'ATGG', 'ATGT', 'ATTA', 'ATTC', 'ATTG', 'ATTT', 'CAAA', 'CAAC', 'CAAG', 'CAAT', 'CACA', 'CACC', 'CACG', 'CACT', 'CAGA', 'CAGC', 'CAGG', 'CAGT', 'CATA', 'CATC', 'CATG', 'CATT', 'CCAA', 'CCAC', 'CCAG', 'CCAT', 'CCCA', 'CCCC', 'CCCG', 'CCCT', 'CCGA', 'CCGC', 'CCGG', 'CCGT', 'CCTA', 'CCTC', 'CCTG', 'CCTT', 'CGAA', 'CGAC', 'CGAG', 'CGAT', 'CGCA', 'CGCC', 'CGCG', 'CGCT', 'CGGA', 'CGGC', 'CGGG', 'CGGT', 'CGTA', 'CGTC', 'CGTG', 'CGTT', 'CTAA', 'CTAC', 'CTAG', 'CTAT', 'CTCA', 'CTCC', 'CTCG', 'CTCT', 'CTGA', 'CTGC', 'CTGG', 'CTGT', 'CTTA', 'CTTC', 'CTTG', 'CTTT', 'GAAA', 'GAAC', 'GAAG', 'GAAT', 'GACA', 'GACC', 'GACG', 'GACT', 'GAGA', 'GAGC', 'GAGG', 'GAGT', 'GATA', 'GATC', 'GATG', 'GATT', 'GCAA', 'GCAC', 'GCAG', 'GCAT', 'GCCA', 'GCCC', 'GCCG', 'GCCT', 'GCGA', 'GCGC', 'GCGG', 'GCGT', 'GCTA', 'GCTC', 'GCTG', 'GCTT', 'GGAA', 'GGAC', 'GGAG', 'GGAT', 'GGCA', 'GGCC', 'GGCG', 'GGCT', 'GGGA', 'GGGC', 'GGGG', 'GGGT', 'GGTA', 'GGTC', 'GGTG', 'GGTT', 'GTAA', 'GTAC', 'GTAG', 'GTAT', 'GTCA', 'GTCC', 'GTCG', 'GTCT', 'GTGA', 'GTGC', 'GTGG', 'GTGT', 'GTTA', 'GTTC', 'GTTG', 'GTTT', 'TAAA', 'TAAC', 'TAAG', 'TAAT', 'TACA', 'TACC', 'TACG', 'TACT', 'TAGA', 'TAGC', 'TAGG', 'TAGT', 'TATA', 'TATC', 'TATG', 'TATT', 'TCAA', 'TCAC', 'TCAG', 'TCAT', 'TCCA', 'TCCC', 'TCCG', 'TCCT', 'TCGA', 'TCGC', 'TCGG', 'TCGT', 'TCTA', 'TCTC', 'TCTG', 'TCTT', 'TGAA', 'TGAC', 'TGAG', 'TGAT', 'TGCA', 'TGCC', 'TGCG', 'TGCT', 'TGGA', 'TGGC', 'TGGG', 'TGGT', 'TGTA', 'TGTC', 'TGTG', 'TGTT', 'TTAA', 'TTAC', 'TTAG', 'TTAT', 'TTCA', 'TTCC', 'TTCG', 'TTCT', 'TTGA', 'TTGC', 'TTGG', 'TTGT', 'TTTA', 'TTTC', 'TTTG', 'TTTT')

fivemers <- c('AAAAA', 'AAAAC', 'AAAAG', 'AAAAT', 'AAACA', 'AAACC', 'AAACG', 'AAACT', 'AAAGA', 'AAAGC', 'AAAGG', 'AAAGT', 'AAATA', 'AAATC', 'AAATG', 'AAATT', 'AACAA', 'AACAC', 'AACAG', 'AACAT', 'AACCA', 'AACCC', 'AACCG', 'AACCT', 'AACGA', 'AACGC', 'AACGG', 'AACGT', 'AACTA', 'AACTC', 'AACTG', 'AACTT', 'AAGAA', 'AAGAC', 'AAGAG', 'AAGAT', 'AAGCA', 'AAGCC', 'AAGCG', 'AAGCT', 'AAGGA', 'AAGGC', 'AAGGG', 'AAGGT', 'AAGTA', 'AAGTC', 'AAGTG', 'AAGTT', 'AATAA', 'AATAC', 'AATAG', 'AATAT', 'AATCA', 'AATCC', 'AATCG', 'AATCT', 'AATGA', 'AATGC', 'AATGG', 'AATGT', 'AATTA', 'AATTC', 'AATTG', 'AATTT', 'ACAAA', 'ACAAC', 'ACAAG', 'ACAAT', 'ACACA', 'ACACC', 'ACACG', 'ACACT', 'ACAGA', 'ACAGC', 'ACAGG', 'ACAGT', 'ACATA', 'ACATC', 'ACATG', 'ACATT', 'ACCAA', 'ACCAC', 'ACCAG', 'ACCAT', 'ACCCA', 'ACCCC', 'ACCCG', 'ACCCT', 'ACCGA', 'ACCGC', 'ACCGG', 'ACCGT', 'ACCTA', 'ACCTC', 'ACCTG', 'ACCTT', 'ACGAA', 'ACGAC', 'ACGAG', 'ACGAT', 'ACGCA', 'ACGCC', 'ACGCG', 'ACGCT', 'ACGGA', 'ACGGC', 'ACGGG', 'ACGGT', 'ACGTA', 'ACGTC', 'ACGTG', 'ACGTT', 'ACTAA', 'ACTAC', 'ACTAG', 'ACTAT', 'ACTCA', 'ACTCC', 'ACTCG', 'ACTCT', 'ACTGA', 'ACTGC', 'ACTGG', 'ACTGT', 'ACTTA', 'ACTTC', 'ACTTG', 'ACTTT', 'AGAAA', 'AGAAC', 'AGAAG', 'AGAAT', 'AGACA', 'AGACC', 'AGACG', 'AGACT', 'AGAGA', 'AGAGC', 'AGAGG', 'AGAGT', 'AGATA', 'AGATC', 'AGATG', 'AGATT', 'AGCAA', 'AGCAC', 'AGCAG', 'AGCAT', 'AGCCA', 'AGCCC', 'AGCCG', 'AGCCT', 'AGCGA', 'AGCGC', 'AGCGG', 'AGCGT', 'AGCTA', 'AGCTC', 'AGCTG', 'AGCTT', 'AGGAA', 'AGGAC', 'AGGAG', 'AGGAT', 'AGGCA', 'AGGCC', 'AGGCG', 'AGGCT', 'AGGGA', 'AGGGC', 'AGGGG', 'AGGGT', 'AGGTA', 'AGGTC', 'AGGTG', 'AGGTT', 'AGTAA', 'AGTAC', 'AGTAG', 'AGTAT', 'AGTCA', 'AGTCC', 'AGTCG', 'AGTCT', 'AGTGA', 'AGTGC', 'AGTGG', 'AGTGT', 'AGTTA', 'AGTTC', 'AGTTG', 'AGTTT', 'ATAAA', 'ATAAC', 'ATAAG', 'ATAAT', 'ATACA', 'ATACC', 'ATACG', 'ATACT', 'ATAGA', 'ATAGC', 'ATAGG', 'ATAGT', 'ATATA', 'ATATC', 'ATATG', 'ATATT', 'ATCAA', 'ATCAC', 'ATCAG', 'ATCAT', 'ATCCA', 'ATCCC', 'ATCCG', 'ATCCT', 'ATCGA', 'ATCGC', 'ATCGG', 'ATCGT', 'ATCTA', 'ATCTC', 'ATCTG', 'ATCTT', 'ATGAA', 'ATGAC', 'ATGAG', 'ATGAT', 'ATGCA', 'ATGCC', 'ATGCG', 'ATGCT', 'ATGGA', 'ATGGC', 'ATGGG', 'ATGGT', 'ATGTA', 'ATGTC', 'ATGTG', 'ATGTT', 'ATTAA', 'ATTAC', 'ATTAG', 'ATTAT', 'ATTCA', 'ATTCC', 'ATTCG', 'ATTCT', 'ATTGA', 'ATTGC', 'ATTGG', 'ATTGT', 'ATTTA', 'ATTTC', 'ATTTG', 'ATTTT', 'CAAAA', 'CAAAC', 'CAAAG', 'CAAAT', 'CAACA', 'CAACC', 'CAACG', 'CAACT', 'CAAGA', 'CAAGC', 'CAAGG', 'CAAGT', 'CAATA', 'CAATC', 'CAATG', 'CAATT', 'CACAA', 'CACAC', 'CACAG', 'CACAT', 'CACCA', 'CACCC', 'CACCG', 'CACCT', 'CACGA', 'CACGC', 'CACGG', 'CACGT', 'CACTA', 'CACTC', 'CACTG', 'CACTT', 'CAGAA', 'CAGAC', 'CAGAG', 'CAGAT', 'CAGCA', 'CAGCC', 'CAGCG', 'CAGCT', 'CAGGA', 'CAGGC', 'CAGGG', 'CAGGT', 'CAGTA', 'CAGTC', 'CAGTG', 'CAGTT', 'CATAA', 'CATAC', 'CATAG', 'CATAT', 'CATCA', 'CATCC', 'CATCG', 'CATCT', 'CATGA', 'CATGC', 'CATGG', 'CATGT', 'CATTA', 'CATTC', 'CATTG', 'CATTT', 'CCAAA', 'CCAAC', 'CCAAG', 'CCAAT', 'CCACA', 'CCACC', 'CCACG', 'CCACT', 'CCAGA', 'CCAGC', 'CCAGG', 'CCAGT', 'CCATA', 'CCATC', 'CCATG', 'CCATT', 'CCCAA', 'CCCAC', 'CCCAG', 'CCCAT', 'CCCCA', 'CCCCC', 'CCCCG', 'CCCCT', 'CCCGA', 'CCCGC', 'CCCGG', 'CCCGT', 'CCCTA', 'CCCTC', 'CCCTG', 'CCCTT', 'CCGAA', 'CCGAC', 'CCGAG', 'CCGAT', 'CCGCA', 'CCGCC', 'CCGCG', 'CCGCT', 'CCGGA', 'CCGGC', 'CCGGG', 'CCGGT', 'CCGTA', 'CCGTC', 'CCGTG', 'CCGTT', 'CCTAA', 'CCTAC', 'CCTAG', 'CCTAT', 'CCTCA', 'CCTCC', 'CCTCG', 'CCTCT', 'CCTGA', 'CCTGC', 'CCTGG', 'CCTGT', 'CCTTA', 'CCTTC', 'CCTTG', 'CCTTT', 'CGAAA', 'CGAAC', 'CGAAG', 'CGAAT', 'CGACA', 'CGACC', 'CGACG', 'CGACT', 'CGAGA', 'CGAGC', 'CGAGG', 'CGAGT', 'CGATA', 'CGATC', 'CGATG', 'CGATT', 'CGCAA', 'CGCAC', 'CGCAG', 'CGCAT', 'CGCCA', 'CGCCC', 'CGCCG', 'CGCCT', 'CGCGA', 'CGCGC', 'CGCGG', 'CGCGT', 'CGCTA', 'CGCTC', 'CGCTG', 'CGCTT', 'CGGAA', 'CGGAC', 'CGGAG', 'CGGAT', 'CGGCA', 'CGGCC', 'CGGCG', 'CGGCT', 'CGGGA', 'CGGGC', 'CGGGG', 'CGGGT', 'CGGTA', 'CGGTC', 'CGGTG', 'CGGTT', 'CGTAA', 'CGTAC', 'CGTAG', 'CGTAT', 'CGTCA', 'CGTCC', 'CGTCG', 'CGTCT', 'CGTGA', 'CGTGC', 'CGTGG', 'CGTGT', 'CGTTA', 'CGTTC', 'CGTTG', 'CGTTT', 'CTAAA', 'CTAAC', 'CTAAG', 'CTAAT', 'CTACA', 'CTACC', 'CTACG', 'CTACT', 'CTAGA', 'CTAGC', 'CTAGG', 'CTAGT', 'CTATA', 'CTATC', 'CTATG', 'CTATT', 'CTCAA', 'CTCAC', 'CTCAG', 'CTCAT', 'CTCCA', 'CTCCC', 'CTCCG', 'CTCCT', 'CTCGA', 'CTCGC', 'CTCGG', 'CTCGT', 'CTCTA', 'CTCTC', 'CTCTG', 'CTCTT', 'CTGAA', 'CTGAC', 'CTGAG', 'CTGAT', 'CTGCA', 'CTGCC', 'CTGCG', 'CTGCT', 'CTGGA', 'CTGGC', 'CTGGG', 'CTGGT', 'CTGTA', 'CTGTC', 'CTGTG', 'CTGTT', 'CTTAA', 'CTTAC', 'CTTAG', 'CTTAT', 'CTTCA', 'CTTCC', 'CTTCG', 'CTTCT', 'CTTGA', 'CTTGC', 'CTTGG', 'CTTGT', 'CTTTA', 'CTTTC', 'CTTTG', 'CTTTT', 'GAAAA', 'GAAAC', 'GAAAG', 'GAAAT', 'GAACA', 'GAACC', 'GAACG', 'GAACT', 'GAAGA', 'GAAGC', 'GAAGG', 'GAAGT', 'GAATA', 'GAATC', 'GAATG', 'GAATT', 'GACAA', 'GACAC', 'GACAG', 'GACAT', 'GACCA', 'GACCC', 'GACCG', 'GACCT', 'GACGA', 'GACGC', 'GACGG', 'GACGT', 'GACTA', 'GACTC', 'GACTG', 'GACTT', 'GAGAA', 'GAGAC', 'GAGAG', 'GAGAT', 'GAGCA', 'GAGCC', 'GAGCG', 'GAGCT', 'GAGGA', 'GAGGC', 'GAGGG', 'GAGGT', 'GAGTA', 'GAGTC', 'GAGTG', 'GAGTT', 'GATAA', 'GATAC', 'GATAG', 'GATAT', 'GATCA', 'GATCC', 'GATCG', 'GATCT', 'GATGA', 'GATGC', 'GATGG', 'GATGT', 'GATTA', 'GATTC', 'GATTG', 'GATTT', 'GCAAA', 'GCAAC', 'GCAAG', 'GCAAT', 'GCACA', 'GCACC', 'GCACG', 'GCACT', 'GCAGA', 'GCAGC', 'GCAGG', 'GCAGT', 'GCATA', 'GCATC', 'GCATG', 'GCATT', 'GCCAA', 'GCCAC', 'GCCAG', 'GCCAT', 'GCCCA', 'GCCCC', 'GCCCG', 'GCCCT', 'GCCGA', 'GCCGC', 'GCCGG', 'GCCGT', 'GCCTA', 'GCCTC', 'GCCTG', 'GCCTT', 'GCGAA', 'GCGAC', 'GCGAG', 'GCGAT', 'GCGCA', 'GCGCC', 'GCGCG', 'GCGCT', 'GCGGA', 'GCGGC', 'GCGGG', 'GCGGT', 'GCGTA', 'GCGTC', 'GCGTG', 'GCGTT', 'GCTAA', 'GCTAC', 'GCTAG', 'GCTAT', 'GCTCA', 'GCTCC', 'GCTCG', 'GCTCT', 'GCTGA', 'GCTGC', 'GCTGG', 'GCTGT', 'GCTTA', 'GCTTC', 'GCTTG', 'GCTTT', 'GGAAA', 'GGAAC', 'GGAAG', 'GGAAT', 'GGACA', 'GGACC', 'GGACG', 'GGACT', 'GGAGA', 'GGAGC', 'GGAGG', 'GGAGT', 'GGATA', 'GGATC', 'GGATG', 'GGATT', 'GGCAA', 'GGCAC', 'GGCAG', 'GGCAT', 'GGCCA', 'GGCCC', 'GGCCG', 'GGCCT', 'GGCGA', 'GGCGC', 'GGCGG', 'GGCGT', 'GGCTA', 'GGCTC', 'GGCTG', 'GGCTT', 'GGGAA', 'GGGAC', 'GGGAG', 'GGGAT', 'GGGCA', 'GGGCC', 'GGGCG', 'GGGCT', 'GGGGA', 'GGGGC', 'GGGGG', 'GGGGT', 'GGGTA', 'GGGTC', 'GGGTG', 'GGGTT', 'GGTAA', 'GGTAC', 'GGTAG', 'GGTAT', 'GGTCA', 'GGTCC', 'GGTCG', 'GGTCT', 'GGTGA', 'GGTGC', 'GGTGG', 'GGTGT', 'GGTTA', 'GGTTC', 'GGTTG', 'GGTTT', 'GTAAA', 'GTAAC', 'GTAAG', 'GTAAT', 'GTACA', 'GTACC', 'GTACG', 'GTACT', 'GTAGA', 'GTAGC', 'GTAGG', 'GTAGT', 'GTATA', 'GTATC', 'GTATG', 'GTATT', 'GTCAA', 'GTCAC', 'GTCAG', 'GTCAT', 'GTCCA', 'GTCCC', 'GTCCG', 'GTCCT', 'GTCGA', 'GTCGC', 'GTCGG', 'GTCGT', 'GTCTA', 'GTCTC', 'GTCTG', 'GTCTT', 'GTGAA', 'GTGAC', 'GTGAG', 'GTGAT', 'GTGCA', 'GTGCC', 'GTGCG', 'GTGCT', 'GTGGA', 'GTGGC', 'GTGGG', 'GTGGT', 'GTGTA', 'GTGTC', 'GTGTG', 'GTGTT', 'GTTAA', 'GTTAC', 'GTTAG', 'GTTAT', 'GTTCA', 'GTTCC', 'GTTCG', 'GTTCT', 'GTTGA', 'GTTGC', 'GTTGG', 'GTTGT', 'GTTTA', 'GTTTC', 'GTTTG', 'GTTTT', 'TAAAA', 'TAAAC', 'TAAAG', 'TAAAT', 'TAACA', 'TAACC', 'TAACG', 'TAACT', 'TAAGA', 'TAAGC', 'TAAGG', 'TAAGT', 'TAATA', 'TAATC', 'TAATG', 'TAATT', 'TACAA', 'TACAC', 'TACAG', 'TACAT', 'TACCA', 'TACCC', 'TACCG', 'TACCT', 'TACGA', 'TACGC', 'TACGG', 'TACGT', 'TACTA', 'TACTC', 'TACTG', 'TACTT', 'TAGAA', 'TAGAC', 'TAGAG', 'TAGAT', 'TAGCA', 'TAGCC', 'TAGCG', 'TAGCT', 'TAGGA', 'TAGGC', 'TAGGG', 'TAGGT', 'TAGTA', 'TAGTC', 'TAGTG', 'TAGTT', 'TATAA', 'TATAC', 'TATAG', 'TATAT', 'TATCA', 'TATCC', 'TATCG', 'TATCT', 'TATGA', 'TATGC', 'TATGG', 'TATGT', 'TATTA', 'TATTC', 'TATTG', 'TATTT', 'TCAAA', 'TCAAC', 'TCAAG', 'TCAAT', 'TCACA', 'TCACC', 'TCACG', 'TCACT', 'TCAGA', 'TCAGC', 'TCAGG', 'TCAGT', 'TCATA', 'TCATC', 'TCATG', 'TCATT', 'TCCAA', 'TCCAC', 'TCCAG', 'TCCAT', 'TCCCA', 'TCCCC', 'TCCCG', 'TCCCT', 'TCCGA', 'TCCGC', 'TCCGG', 'TCCGT', 'TCCTA', 'TCCTC', 'TCCTG', 'TCCTT', 'TCGAA', 'TCGAC', 'TCGAG', 'TCGAT', 'TCGCA', 'TCGCC', 'TCGCG', 'TCGCT', 'TCGGA', 'TCGGC', 'TCGGG', 'TCGGT', 'TCGTA', 'TCGTC', 'TCGTG', 'TCGTT', 'TCTAA', 'TCTAC', 'TCTAG', 'TCTAT', 'TCTCA', 'TCTCC', 'TCTCG', 'TCTCT', 'TCTGA', 'TCTGC', 'TCTGG', 'TCTGT', 'TCTTA', 'TCTTC', 'TCTTG', 'TCTTT', 'TGAAA', 'TGAAC', 'TGAAG', 'TGAAT', 'TGACA', 'TGACC', 'TGACG', 'TGACT', 'TGAGA', 'TGAGC', 'TGAGG', 'TGAGT', 'TGATA', 'TGATC', 'TGATG', 'TGATT', 'TGCAA', 'TGCAC', 'TGCAG', 'TGCAT', 'TGCCA', 'TGCCC', 'TGCCG', 'TGCCT', 'TGCGA', 'TGCGC', 'TGCGG', 'TGCGT', 'TGCTA', 'TGCTC', 'TGCTG', 'TGCTT', 'TGGAA', 'TGGAC', 'TGGAG', 'TGGAT', 'TGGCA', 'TGGCC', 'TGGCG', 'TGGCT', 'TGGGA', 'TGGGC', 'TGGGG', 'TGGGT', 'TGGTA', 'TGGTC', 'TGGTG', 'TGGTT', 'TGTAA', 'TGTAC', 'TGTAG', 'TGTAT', 'TGTCA', 'TGTCC', 'TGTCG', 'TGTCT', 'TGTGA', 'TGTGC', 'TGTGG', 'TGTGT', 'TGTTA', 'TGTTC', 'TGTTG', 'TGTTT', 'TTAAA', 'TTAAC', 'TTAAG', 'TTAAT', 'TTACA', 'TTACC', 'TTACG', 'TTACT', 'TTAGA', 'TTAGC', 'TTAGG', 'TTAGT', 'TTATA', 'TTATC', 'TTATG', 'TTATT', 'TTCAA', 'TTCAC', 'TTCAG', 'TTCAT', 'TTCCA', 'TTCCC', 'TTCCG', 'TTCCT', 'TTCGA', 'TTCGC', 'TTCGG', 'TTCGT', 'TTCTA', 'TTCTC', 'TTCTG', 'TTCTT', 'TTGAA', 'TTGAC', 'TTGAG', 'TTGAT', 'TTGCA', 'TTGCC', 'TTGCG', 'TTGCT', 'TTGGA', 'TTGGC', 'TTGGG', 'TTGGT', 'TTGTA', 'TTGTC', 'TTGTG', 'TTGTT', 'TTTAA', 'TTTAC', 'TTTAG', 'TTTAT', 'TTTCA', 'TTTCC', 'TTTCG', 'TTTCT', 'TTTGA', 'TTTGC', 'TTTGG', 'TTTGT', 'TTTTA', 'TTTTC', 'TTTTG', 'TTTTT')


read.model.file <- function(model.file, variantcolumn=FALSE){
  data <- read.table(model.file)
  if(variantcolumn){
    colnames(data) <- c("kmer",  "variant",	"level_mean",	"level_stdv",	"sd_mean",	"sd_stdv",	"weight")
  }else{
    colnames(data) <- c("kmer", "level_mean",	"level_stdv",	"sd_mean", "sd_stdv",	"weight")
  }
  return(data)
}

read.events.file <- function(events.file, input.events=FALSE, numRows=NA){
  #numRows -- max amount of rows to read in -- shoot higher - speeds up
  ## if numRows given -- it assumes large fle and also attempts to use colclasses
  if(is.na(numRows)){data <- read.table(events.file)}
  else {data <- read.table(events.file, nrows = numRows, colClasses = c("numeric", "numeric", "numeric", "numeric"))}     
  if(input.events){
    colnames(data) <- c("mean", "stddev",  "start", "length")
  } else {##template or comp events
    colnames(data) <- c("mean", "stddev",  "start", "length", "model_state",	"model_level",	"move",	"p_model_state", "mp_state",	"p_mp_state",	"p_A",	"p_C",	"p_G",	"p_T")
  }
  return(data)
}


####### ALIGN
library(Biostrings)
# nw.align <- function(s1, s2, d){
#   # s1 = sequence 1
#   # s2 = sequence 2
#   # d = gap penalty
#   n = nchar(s1)
#   m = nchar(s2)
#   A = matrix(data = rep(0, (n+1)*(m+1)), nrow = n+1)
#   P = matrix(data = rep(0, (n+1)*(m+1)), nrow = n+1)
#   diagonal <- 0
#   up <- 1
#   left <- 2
#   ## init -- fill in left and top
#   for(i in 2:(n+1)){A[i,0] <- (i-1)*d; P[i,0] <- left}
#   for(j in 2:(m+1)){A[0,j] <- (j-1)*d}
#   ## init pointer
#   
#   ## iter from top right to bottom left
#   for(i in 2:(n+1)){
#     for(j in 2:(m+1)){
#       match = -1 + (substr(s1,i,i) == substr(s2,j,j))
#     }
#   }
# }
phmm.align <- function(s1, s2, logprobs=TRUE){
  ## logprobs = whether or not transition and initial matrices are logged yet, assumes false
  if(!(logprobs)){
    initial = log(initial)
    transitions = log(transitions)
  }
  num.states <- length(states)
  num.emits <- length(emitted.data)
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)   
  ## need to add log_probs instead of multiply probs to prevent underflow
  
  Viterbi[ ,1] <- initial + dnorm(emitted.data[1], mean = emissions[1, ], sd = emissions[2, ], log = TRUE)
  pointer[1, ] <- 1
  
  for (j in 2:num.emits){
    #     selection <- Viterbi[,j-1] + transitions
    #     maxstates <- apply(selection, 2, which.max)
    #     maxes <- apply(selection, 2, max)
    #     Viterbi[ ,j] <- dnorm(emitted.data[j], mean = emissions[1, ], sd = emissions[2, ], log = TRUE) + maxes
    #     pointer[j, ] <- maxstates
    for (i in 1:num.states){
      selection <- Viterbi[,j-1] + transitions[,i]
      maxstate <- which.max(selection)
      Viterbi[i,j] <- dnorm(emitted.data[j], mean = emissions[1, i], sd = emissions[2, i], log = TRUE) + selection[maxstate]
      pointer[j,i] <- maxstate 
    }
  }     
  
  viterbi_path <- matrix(data = rep(0, num.emits), nrow = 1)
  viterbi_path[num.emits] <- which.max(Viterbi[,num.emits]);  
  viterbi_prob <- Viterbi[viterbi_path[num.emits], num.emits]
  
  for (j in num.emits:2){
    viterbi_path[j-1] <- pointer[j,viterbi_path[j]]
  }
  seqinfo=getseqfxn(states, viterbi_path)
  return(list(viterbi_path=viterbi_path, viterbi_seq=seqinfo$seq, viterbi_prob=viterbi_prob, moves=seqinfo$moves))
}




############      Dec 15, 2015      ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
############ for DNA puff puff pass ############
viterbi.puff <- function(emissions, transitions, initial, states, emitted.data, emodel="normal", getseqfxn=get.sequence, logprobs=FALSE){
  ## logprobs = whether or not transition and initial matrices are logged yet, assumes false
  if(!(logprobs)){
    initial = log(initial)
    transitions = log(transitions)
  }
  num.states <- length(states)
  num.emits <- length(emitted.data) 
  ## need to add log_probs instead of multiply probs to prevent underflow
  
  if (emodel == "normal"){V=viterbi.normal(emissions, transitions, initial, emitted.data, num.states, num.emits)}
  else if (emodel == "poisson"){V=viterbi.poisson(emissions, transitions, initial, emitted.data, num.states, num.emits)}
  else if (emodel == "exponential"){V=viterbi.exponential(emissions, transitions, initial, emitted.data, num.states, num.emits)}
  else if (emodel == "geometric"){V=viterbi.geometric(emissions, transitions, initial, emitted.data, num.states, num.emits)}  
  
  viterbi_path <- matrix(data = rep(0, num.emits), nrow = 1)
  viterbi_path[num.emits] <- which.max(V$Viterbi[,num.emits]);  
  viterbi_prob <- V$Viterbi[viterbi_path[num.emits], num.emits]
  
  for (j in num.emits:2){
    viterbi_path[j-1] <- V$pointer[j,viterbi_path[j]]
  }
  seqinfo=getseqfxn(states, viterbi_path)
  return(list(viterbi_path=viterbi_path, viterbi_seq=seqinfo$seq, viterbi_prob=viterbi_prob, moves=seqinfo$moves))
}

##NORMAL
viterbi.normal <- function(emissions, transitions, initial, emitted.data, num.states, num.emits){
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)  
  Viterbi[ ,1] <- initial + dnorm(emitted.data[1], mean = emissions[1, ], sd = emissions[2, ], log = TRUE)
  pointer[1, ] <- 1
  f <- function(x){i <- which.max(x); y <- x[i]; return(c(i,y))}
  for (j in 2:num.emits){
    selection <- Viterbi[,j-1] + transitions
    for (i in 1:num.states){
      maxstate <- which.max(selection[,i])
      Viterbi[i,j] <- dnorm(emitted.data[j], mean = emissions[1, i], sd = emissions[2, i], log = TRUE) + selection[maxstate,i]
      pointer[j,i] <- maxstate 
    }
  }   
  return(list(Viterbi=Viterbi, pointer=pointer))
}

##POISSON
viterbi.poisson <- function(emissions, transitions, initial, emitted.data, num.states, num.emits){
  ## rounds floats to integers
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)  
  Viterbi[ ,1] <- initial + dpois(round(emitted.data[1]), lambda = emissions[1, ], log=TRUE)
  pointer[1, ] <- 1
  f <- function(x){i <- which.max(x); y <- x[i]; return(c(i,y))}
  for (j in 2:num.emits){
    selection <- Viterbi[,j-1] + transitions
    for (i in 1:num.states){
      maxstate <- which.max(selection[,i])
      Viterbi[i,j] <- dpois(round(emitted.data[j]), lambda = emissions[1, i], log = TRUE) + selection[maxstate,i]
      pointer[j,i] <- maxstate 
    }
  }  
  return(list(Viterbi=Viterbi, pointer=pointer))
}

viterbi.exponential <- function(emissions, transitions, initial, emitted.data, num.states, num.emits){
  ## rounds floats to integers
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)  
  Viterbi[ ,1] <- initial + dexp(emitted.data[1], rate = emissions[1, ], log=TRUE)
  pointer[1, ] <- 1
  f <- function(x){i <- which.max(x); y <- x[i]; return(c(i,y))}
  for (j in 2:num.emits){
    selection <- Viterbi[,j-1] + transitions
    for (i in 1:num.states){
      maxstate <- which.max(selection[,i])
      Viterbi[i,j] <- dexp(emitted.data[j], rate = emissions[1, i], log = TRUE) + selection[maxstate,i]
      pointer[j,i] <- maxstate 
    }
  }  
  return(list(Viterbi=Viterbi, pointer=pointer))
}

viterbi.geometric <- function(emissions, transitions, initial, emitted.data, num.states, num.emits){
  ## rounds floats to integers
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)  
  Viterbi[ ,1] <- initial + dgeom(round(emitted.data[1]), prob = emissions[1, ], log=TRUE)
  pointer[1, ] <- 1
  f <- function(x){i <- which.max(x); y <- x[i]; return(c(i,y))}
  for (j in 2:num.emits){
    selection <- Viterbi[,j-1] + transitions
    for (i in 1:num.states){
      maxstate <- which.max(selection[,i])
      Viterbi[i,j] <- dgeom(round(emitted.data[j]), prob = emissions[1, i], log = TRUE) + selection[maxstate,i]
      pointer[j,i] <- maxstate 
    }
  }  
  return(list(Viterbi=Viterbi, pointer=pointer))
}

#### VECTORIZED FORWARDS ########
forward.puff <- function(emissions, transitions, initial, states, emitted.data, emodel="normal"){
  num.states <- length(states)
  num.emits <- length(emitted.data)
  Forward <- matrix(data = rep(0, num.states*num.emits), nrow = num.states)
  scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
  
  #model
  if (emodel == "normal"){emodel.fxn <- puff.normal}
  else if (emodel == "exponential"){emodel.fxn <- puff.exponential}
  else if (emodel == "poisson"){emodel.fxn <- puff.poisson}
  else if (emodel == "geometric"){emodel.fxn <- puff.geometric}
  
  ## initial
  Forward[, 1] <- initial*emodel.fxn(emitted.data[1], emissions)
  ## scale to prevent underflow -- keep track of scaling
  scalefactors[1,1] <- sum(Forward[, 1])
  scalefactors[2,1] <- log(scalefactors[1,1])
  Forward[,1] <- Forward[,1]/scalefactors[1,1]
  
  ## iterate
  for(k in 2:num.emits){
    emit <- emodel.fxn(emitted.data[k], emissions)
    Forward[, k] <- emit* Forward[,k-1] %*% transitions ## same as emit* Forward[,k-1] * colSums(transitions)
    scalefactors[1,k] <- sum(Forward[, k])
    scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k-1]
    Forward[,k] <- Forward[,k]/scalefactors[1,k]
  }
  
  return(list(forward=Forward, scales=scalefactors))
  ## mutiply forward column by row2,samecol in scale factors OR by exp(row3,samecol in scalfactors) to get actual value for forward
  ## 
}




#### VECTORIZED BACKWARDS ########
backward.puff <- function(emissions, transitions, initial, states, emitted.data, emodel="normal"){
  num.states <- length(states)
  num.emits <- length(emitted.data)
  Backward = matrix(data = rep(0, num.states*num.emits), nrow = num.states)
  scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
  
  #model
  if (emodel == "normal"){emodel.fxn <- puff.normal}
  else if (emodel == "exponential"){emodel.fxn <- puff.exponential}
  else if (emodel == "poisson"){emodel.fxn <- puff.poisson}
  else if (emodel == "geometric"){emodel.fxn <- puff.geometric}
  
  ## initial
  Backward[ , num.emits] <- 1
  
  ## scale to prevent underflow -- keep track of scaling
  scalefactors[1,num.emits] <- sum(Backward[, num.emits])
  scalefactors[2,num.emits] <- log(scalefactors[1,num.emits])
  Backward[,num.emits] <- Backward[,num.emits]/scalefactors[1,num.emits]
  
  ## iterate
  for(k in (num.emits-1):1){
#     emit <- matrix(dnorm(emitted.data[k+1], mean = emissions[1, ], sd = emissions[2, ]))
    emit <- matrix(emodel.fxn(emitted.data[k+1], emissions))
    #     print(Backward[, k+1] * emit)
    Backward [, k] <- transitions %*% (Backward[, k+1] * emit)
    scalefactors[1,k] <- sum(Backward[, k])
    scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k+1]
    Backward[,k] <- Backward[,k]/scalefactors[1,k]
  }
  return(list(backward=Backward, scales=scalefactors))
}


puff.normal <- function(x, emissions){
  dnorm(x, mean = emissions[1, ], sd = emissions[2, ], log=FALSE)
}

puff.exponential <- function(x, emissions){
  dexp(x = x, rate = emissions[1, ], log = FALSE)
}

puff.poisson <- function(x, emissions){
  dpois(x = round(x), lambda = emissions[1, ], log = FALSE)
}

puff.geometric <- function(x, emissions){
  dgeom(x = round(x), prob = emissions[1, ], log = FALSE)
}
