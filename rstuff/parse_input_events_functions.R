
find_position_in_input <- function(subset, input, window=50){
  ##subset can be wither template or complement events from ONT basecaller
  currentpos <- 0
  diff <- 1
  in.size <- length(input)
  sub.size <- length(subset)
  lastpos <- in.size - window + 1
  subset.win <- subset[1:window]
#   print(subset.win)
#   print(subset.win - input[51:(51+50-1)])
  while( diff != 0 && currentpos <= lastpos ){
    currentpos <- currentpos + 1
#     print(currentpos)
    diff <- sum(subset.win - input[currentpos:(currentpos+window-1)])
#     print(diff)
  }
# print(diff)
#   print(c(diff,currentpos))
  if(diff == 0){
    if(sum(subset-input[currentpos:(currentpos+sub.size-1)]) == 0){
      ALL=TRUE
    } else {ALL=FALSE}
  } else {
    return("Did not find positions")
  }
  return(list(start=currentpos, end=currentpos+sub.size-1, all=ALL))
}


max_window_index <-function(vec, winsize=60){
  endpos <- length(vec)-winsize+1
  score <- sum(vec[1:winsize])
  maxsum <- score
  maxindex <- 1
  for(i in 2:endpos){
    score <- score - (vec[i-1]) + vec[i+winsize-1]
    if(score > maxsum){
      maxsum <- score
      maxindex <- i
    }
  }
  return(list(max.score=maxsum, max.index=maxindex))
}


process_events <- function(t_events, c_events, i_events){
  template <- read.events.file(t_events, input.events = FALSE)
  complement <- read.events.file(c_events, input.events = FALSE)
  input <- read.events.file(i_events, input.events = TRUE)
  t_on_i <- find_position_in_input(subset=template$mean, input=input$mean, window=50) # 51, 3473, FALSE
  end.template <- template$mean[(length(template$mean)-20):length(template$mean)] ## 
  ## get actual end of template in input
  tend_on_i <- find_position_in_input(subset=end.template, input=input$mean, window=10) # 3463, 3483, TRUE
  # template = 51 - 3483 -- with 10 events erased from input (probably 10 0 moves)
  c_on_i <- find_position_in_input(subset=complement$mean, input=input$mean, window=50) # 3544, 6260, FALSE
  end.comp <- complement$mean[(length(complement$mean)-20):length(complement$mean)] 
  cend_on_i <- find_position_in_input(subset=end.comp, input=input$mean, window=5) # starts at 6247 and ends on 6267 in input
  ## hairpin length = (actual beginning of complement -1) - (actual end of template + 1) + 1 = actBeginComp - ActEndTemp - 1
  hpstart <- tend_on_i$end + 1
  ###TODO
  hplen <- c_on_i$begin - tend_on_i$end - 1 
  # is max signal out of all input events inside hairpin?
  max_in_hp = which.max(input$mean) %in% (tend_on_i$end+1):(c_on_i$begin-1)
  rel_index_of_max_in_hp <- which.max(input$mean[(tend_on_i$end+1):(c_on_i$begin-1)]) 
  return(list(hpstart=hpstart, hpend=hpend, max_in_hp=max_in_hp, ))
}


lead_hmm <- function(L, data){
#   L=list(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14)
  ## state 15 is end state
  emit = matrix(rep(0, 15*2), nrow=2)
  emit[1,] = c(sapply(L, mean), 60)
  emit[2,] = c(sapply(L, sd), 5)
  init = matrix(c(0.4,0.3,0.2,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)+0.001
  init = init/sum(init)
  t = matrix(rep(0, 15*15), nrow = 15)
  t[15,15] <- 1.0
  for(i in 1:14){t[i,i] <- 0.3}
  for(i in 1:13){t[i,i+1] = 0.35}
  for(i in 1:12){t[i,i+2] = 0.2}
  for(i in 1:11){t[i,i+3] = 0.1}
  for(i in 1:10){t[i,i+4] = 0.001}
  for(i in 1:9){t[i,i+5] = 0.001}
  for(i in 1:8){t[i,i+6] = 0.001}
  for(i in 1:7){t[i,i+7] = 0.001}
  t[12,15] <- 0.05
  t[13,15] <- 0.1
  t[14,15] <- 0.2 ## this was 0.35 before this line -- should I really make this change?
  for(i in 1:14){t[i,] <- t[i,]/sum(t[i,])}
  v<-viterbi(emissions = emit, transitions = t, initial = init, states = 1:15, emitted.data = data)
  print(v$viterbi_path)
}


hp_hmm <- function(L, data){
  ## state B,1,2,3,4,5,E = 7 states
  emit = matrix(rep(0, 7*2), nrow=2)
  emit[1,] = c(65, sapply(L, mean), 65)
  emit[2,] = c(6, sapply(L, sd), 6)
  init = matrix(c(0.7,0.2,0.1,0,0,0,0),nrow=1)
  init = init/sum(init)
  t = matrix(rep(0, 7*7), nrow = 7)
  t[7,7] <- 1.0
  for(i in 1:7){t[i,i] <- 0.3}
  for(i in 1:6){t[i,i+1] = 0.35}
  for(i in 1:5){t[i,i+2] = 0.2}
  for(i in 1:4){t[i,i+3] = 0.1}
  for(i in 1:3){t[i,i+4] = 0.001}
  for(i in 1:2){t[i,i+5] = 0.001}
  t[4,7] <- 0.05
  t[5,7] <- 0.1
  t[5,6] <- 0.7 ## state 5 usually goes directly to 6 (occasionally 2 events, but have not seen more)
  for(i in 1:7){t[i,] <- t[i,]/sum(t[i,])}
  v<-viterbi(emissions = emit, transitions = t, initial = init, states = 1:7, emitted.data = data)
  return(list(statepath=v$viterbi_path,logvitprob=v$viterbi_prob))
}

hp.g2=c(100.29852151499328,99.161932576497392, 102.85814409304886, 97.273369750976556, 84.540624389648428, 94.787215795023016, 89.73405826344208, 91.323413899739577, 84.100183512369782)
hp.g3=c(112.33164294433594, 114.80032583797677, 117.82765869140624, 107.18157348632812, 110.93718450881141, 105.17127548217773, 117.23543159179687, 124.89072875976562, 119.97993881835937, 109.34451237605168, 129.74444488525393, 129.65196413730055, 131.34582204718339, 121.17269575517254, 123.47322377472675, 129.76626403808595, 110.33153889973957, 119.77973129734848, 120.59478560353053, 125.10305627893518, 109.81249595528217, 108.78670242309569, 103.07503622731855, 104.54422299985531, 123.57926435786523, 128.20344676050647, 113.73188480377196)
hp.g4=c(92.982586669921872, 103.44664348098466, 104.77688069661457, 99.539638264973945, 104.21992318725586, 102.61014038085936, 97.144240493774419)
hp.g5=c(56.96401758829753,63.484381245457847, 62.753639346852012, 58.064201733398434, 63.396030273437496, 61.145037554572603, 66.28625891397165, 67.992860229492194, 52.73078303019205, 76.449917918238143, 55.214107560797352, 55.413723179659392, 53.806870727539057, 57.419842170266548, 57.461720956089984,  54.584843515249396)
hp.g6=c(47.020195896474561, 47.536156372923955, 46.816061323280429, 46.907960586547844, 45.953837832496284, 46.560520554001272, 45.435220301011029, 46.938543456346423, 46.468505504052963, 46.730685253425186, 47.970399676067068, 47.065632857771114, 46.620575853837401, 45.442259023813101, 48.872249058314729, 46.97980486012468, 45.1288427474651, 46.615050863125525, 48.514198106695851, 44.32685353246228, 43.606539378086453, 43.78506983263739, 43.807955709706178, 43.65661553440323, 45.06628741096047, 49.930329463930633, 51.852153757579288)
hp.L=list(hp.g2,hp.g3,hp.g4,hp.g5,hp.g6)

hp_hmm_loop <- function(data, winsize=60, step=15, L){
  ## data vector of means
  ## winsize -- e,g, 60
  ## step -- moving sliding window
  ## L is a list of groups where mu and sigma can be found for profile hmm at each step -- e.g. hp.L
  start <- 1
  end <- length(data)-60+1
  ans <- list(statepath=1, logvitprob=-1e9)
  winstart <- 1
  for(i in seq(start, end, step)){
    vec <- data[i:(i+winsize-1)]
    currans <- hp_hmm(L=L, data=vec)
    if(currans$logvitprob > ans$logvitprob){
      ans <- currans
      winstart <- i
      print(list(currbest=ans, winstart=winstart))
    }
  }
  return(list(statepath=ans$statepath, logvitprob=ans$logvitprob, winstart=winstart))
}
