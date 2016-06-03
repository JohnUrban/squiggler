## set pwd to ~/path/to/squiggler/example_dtw/ 
setwd("~/searchPaths/bitbucket/squiggler/example_dtw/") 
source("~/searchPaths/bitbucket/squiggler/example_dtw/hmm_functions.R")

## toy example

## emissions: mean and sd
emissions = as.matrix(data.frame(A=c(10,sqrt(5)), C=c(40,sqrt(5)), G=c(70,sqrt(5)), T=c(100,sqrt(5))))
# emissions
emissions = as.matrix(data.frame(A=c(40,2), C=c(50,2.5), G=c(60,2.2), T=c(70,1.8)))

##transitions
a = c(0.1,  0.2,	0.3,	0.4)
c = c(0.4,  0.3,	0.2,	0.1)
g = c(0.25,  0.25,	0.15,	0.35)
t = c(0.3,  0.2,	0.3,	0.2)
uniform = c(0.25,0.25,0.25,0.25)
transitions = matrix(data = c(a,c,g,t), nrow = 4, byrow = TRUE)
transitions = matrix(data = rep(uniform, 4), nrow = 4, byrow = TRUE)
# transitions

initial = c(0.25,0.25,0.25,0.25)
states = c("A","C","G","T")

## examples of generating markov chain or hidden markov chain -- short 10-200 nt
# generate.statepath(transitions, initial, states, length=10)
# generate.emissions(emissions, transitions, initial, states, length=10)

answer <- generate.emissions(emissions, transitions, initial, states, length=10) 
# answer
v <- viterbi(emissions, transitions, initial, states, answer$emitted.data)
# v
compare.statepath(v$viterbi_path, answer$statepath)

f <- forward(emissions, transitions, initial, states, answer$emitted.data)
b <- backward(emissions, transitions, initial, states, answer$emitted.data)
d <- prob.data(f)
sum(f[,1]*b[,1])/d
c <- centroid(f,b,states)
compare.statepath(c$centroid.path, answer$statepath)


v$viterbi_seq
c$centroid.seq
answer$statepathseq

### Example with longer sequences > 200 nt -- F, B and V algorithms need log and scaling procedures to prevent underflow
### I only wrote the log version of the viterbi
source("~/searchPaths/bitbucket/squiggler/example_dtw/hmm_functions.R")
answer <- generate.emissions(emissions, transitions, initial, states, length=10) 
fl <- forward.long(emissions, transitions, initial, states, answer$emitted.data)
dl <- prob.data.long(Forward = fl$forward, scalefactors = fl$scales)
bl <- backward.long(emissions, transitions, initial, states, answer$emitted.data)
v <- viterbi(emissions, transitions, initial, states, answer$emitted.data)
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states)
compare.statepath(v$viterbi_path, answer$statepath)
compare.statepath(cl$centroid.path, answer$statepath)
compare.seq(v$viterbi_seq, answer$statepathseq)
compare.seq(cl$centroid.seq, answer$statepathseq)
####################################################################################
## more real example
emissions = data.frame(A=c(60,sqrt(5)), C=c(65,sqrt(10)), G=c(70,sqrt(20)), T=c(55,5))
a = c(0.1,  0.2,	0.3,	0.4)
c = c(0.4,  0.3,	0.2,	0.1)
g = c(0.25,  0.25,	0.15,	0.35)
t = c(0.3,  0.2,	0.3,	0.2)
transitions = matrix(data = c(a,c,g,t), nrow = 4, byrow = TRUE)
initial = c(0.25,0.25,0.25,0.25)
states = c("A","C","G","T")
answer <- generate.emissions(emissions, transitions, initial, states, length=5000) 
fl <- forward.long(emissions, transitions, initial, states, answer$emitted.data)
dl <- prob.data.long(Forward = fl$forward, scalefactors = fl$scales)
bl <- backward.long(emissions, transitions, initial, states, answer$emitted.data)
v <- viterbi(emissions, transitions, initial, states, answer$emitted.data)
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states)
compare.statepath(v$viterbi_path, answer$statepath)
compare.statepath(cl$centroid.path, answer$statepath) 
compare.seq(v$viterbi_seq, answer$statepathseq)
compare.seq(cl$centroid.seq, answer$statepathseq)

####################################################################################
## sampling
simulateACGT(n=100, seqlen=100)
## 1/16/15: in 100 trials generating seqs of length 100 with diff params each time, centroid got 
##  a higher pct identity to the true sequence than viterbi 62% of the time (Viterbi won 38%).
## The centroid might be an improvement to the basecaller which I believe uses viterbi




###########
source("~/searchPaths/bitbucket/squiggler/example_dtw/hmm_functions.R")
## made function to automatically generate transition probs
## put dimers, trimers, fourmers, and fivemers in environment
#Dimers:
dimer.t <- generate.kmer.transition.probs(dimers)
dimer.e <- generate.kmer.emission.probs(dimers)
dimer.i <- generate.kmer.initial.probs(dimers)
answer <- generate.emissions(emissions = dimer.e, transitions = dimer.t, initial = dimer.i, states = dimers, length=10000) 


fl <- forward.long(emissions = dimer.e, transitions = dimer.t, initial = dimer.i, states = dimers, answer$emitted.data)
dl <- prob.data.long(Forward = fl$forward, scalefactors = fl$scales)
bl <- backward.long(emissions = dimer.e, transitions = dimer.t, initial = dimer.i, states = dimers, answer$emitted.data)
v <- viterbi(emissions = dimer.e, transitions = dimer.t, initial = dimer.i, states = dimers, answer$emitted.data)
# v.f <- viterbi.faster(emissions = dimer.e, transitions = dimer.t, initial = dimer.i, states = dimers, answer$emitted.data)
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states = dimers, getseqfxn = get.sequence)
compare.statepath(v$viterbi_path, answer$statepath)
compare.statepath(cl$centroid.path, answer$statepath)   
compare.statepath(cl$centroid.path, v$viterbi_path)
#compare seqs
compare.seq(v$viterbi_seq, answer$statepathseq)
compare.seq(cl$centroid.seq, answer$statepathseq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)

nchar(v$viterbi_seq) == nchar(cl$centroid.seq)
v$viterbi_seq
cl$centroid.seq


## FIVEMERS ###############################################
source("~/searchPaths/bitbucket/squiggler/example_dtw/hmm_functions.R")
fivemer.t <- generate.kmer.transition.probs(fivemers)
fivemer.e <- generate.kmer.emission.probs(fivemers)
fivemer.i <- generate.kmer.initial.probs(fivemers)
answer <- generate.emissions(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, length=500) 

date()
fl <- forward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, answer$emitted.data)
date()
dl <- prob.data.long(Forward = fl$forward, scalefactors = fl$scales)
date()
bl <- backward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, answer$emitted.data)
date()
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, answer$emitted.data)
# v.f <- viterbi.faster(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, answer$emitted.data)
date()
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states = fivemers)
date()
compare.statepath(v$viterbi_path, answer$statepath)
# compare.statepath(v.f$viterbi_path, answer$statepath)
compare.statepath(cl$centroid.path, answer$statepath)  
compare.statepath(cl$centroid.path, v$viterbi_path)
# compare.statepath(v$viterbi_path, v.f$viterbi_path)
#compare seqs
compare.seq(v$viterbi_seq, answer$statepathseq)
compare.seq(cl$centroid.seq, answer$statepathseq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)

## 10 bp seq -- 
## F: ~144 MB memory, 98-100% CPU, 1 minute 52 sec (112 seconds -- 11.2 seconds per base (196 hours for 63k read))
## prob.data -- instantaneous
## B: ~144 MB memory, 98-100% CPU, 1 minute, 55 sec-- roughly same as forward
## V: ~1 seconds
## C: < 1 second


## HMM1 -- real data -- no bias correction -- matches what ONT puts out 50-75%
## NOTE: agreement with ONT changes each time b/c transition probs change each time.
## I can learn their transition probabilities by maximizing an alignment
## SCIARA1 - HMM1 - 54.7% ID with ONT
fivemer.t <- generate.kmer.transition.probs(fivemers)
model_r7.3 <- read.model.file("template_model_r7.3.tsv")
fivemer.e <- t(model_r7.3[,2:3])
data <- read.events.file("../data/events/sciara_ch111_file162_template.tsv")
emitted.data <- data$mean
sci_ch111_template="CTGGGTGTTTATACCAACGTGTTCTACAACGCATTATGCAGTATAATCCTCCATTTCTTAGTCGCCTTTATTACTAGAAATGTTTTCTAGTGAGGGATTCTTTCGGGCCCATGCCGAAATGAAACCAGCATTTCAAAGCGCTACCCCGCTCATGCTCGTATTTGTCCCGTGCATA"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
fl <- forward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
bl <- backward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states = fivemers)


compare.seq(substr(sci_ch111_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(sci_ch111_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)
compare.seq(v$viterbi_seq, cl$centroid.seq)

##SCIARA2 - HMM1 - 72.2% ID with ONT
data <- read.events.file("../data/events/sciara_ch224_file20_template.tsv")
emitted.data <- data$mean
sci_ch224_template="AGTATAGTCTCATCGGGTTAAGAAAATATTCTCGCAGGTTAGTTGTGTTTCTGTTAATACCATTTTACTACGACTTCGAGCCGTGAGGCCCGGGTGAGCTGGGCACAAGGATAGGAGCGACGTTCCTCTTCCACCTTAATCCTCCAGTGGTAGCCATTAGCCATGATGTAACAAT"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
fl <- forward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
bl <- backward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states = fivemers)
compare.seq(substr(sci_ch224_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(sci_ch224_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)
compare.seq(v$viterbi_seq, cl$centroid.seq)


##LAMBDA - HMM1 - 60.7% ID with ONT
model_r7 <- read.model.file("template_model_r7.tsv")
fivemer.e <- t(model_r7[,2:3])
fivemer.e <- t(model_r7.3[,2:3]) ## r7.3 actually works better
data <- read.table("../data/events/lambda_ch230_file5_template.tsv")
emitted.data <- data$V1 
lambda_ch230_template="TTTTCAATTTATGAAGGTCCGTCAGGCGTTCGTTCTTCTTCGTCATATCCTGGTTTCTTCATGAGCTCAGGGAAGGAGCACAATCATGAAGCATTTTGACATATGTCGTTCCTGAGTTTTGTCCTGTCGCATCATAATGGAAATCCTAGGGCAAACGGGCTTTCGGTGCCAATAT"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
fl <- forward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
bl <- backward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states = fivemers)
compare.seq(substr(lambda_ch230_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(lambda_ch230_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)
compare.seq(v$viterbi_seq, cl$centroid.seq)
v$moves
data$V7[1:50]

#LAMBDA2 - HMM1 - 53.7% ID with ONT
data <- read.events.file("../data/events/lambda_ch304_file2_template.tsv")
emitted.data <- data$mean
lambda_ch304_template <- "GTGTTATGCTGATCCCATTATGGAAGTTTCATTAGAGACAAAGCCGTCGGGACTATCATTATCGATACTGAACTCCCAACGTGCGTGGAGCATCAACCGTTACGGCTTAATTATGAAAATACTACGTATTAATTATTATCCATCTATCTAACGTAAAATAGGTGCGACATCAAGT"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50], getseqfxn = get.sequence.withgaps)
compare.seq(substr(lambda_ch304_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
v$moves
data$move[1:50]

### Note in squiggler, we should make sure that calculating the F, B, and V are parallel-izable.
### For this 62-63k read, running this on the GUI, my machine would take 196 hours for forward and backward.. too much
### Therefore, the centroid is hard to compute
### seemed CPU intensive, but not bad on memory ~ 600-700 MB for the Forward
## viterbi can be computer quicker -- went up 1 > 1.15 GB mem:

## the seq is off -- mess around with scaling
## scaling:
##ONT says:
# In the read files there are some attributes in UniqueGlobalKey/channel_id. The key ones are digitisation (in adc units), offset (in adc units) and range (in pA). That gives the transformation from the raw adc units to pA.
# pA = (raw + offset)*range/digitisation
# These attributes are also present in the bulk fast5, associated with the raw trace in Raw/Channel_X/Meta

sci_230_offset = 5
sci_230_range = 1861.82
sci_230_digitization = 8192
pA = (emitted.data + sci_230_offset)*sci_230_range/sci_230_digitization
head(pA)
template="TTTTCAATTTATGAAGGTCCGTCAGGCGTTCGTTCTTCTTCGTCATATCCTGGTTTCTTCATGAGCTCAGGGAAGGAGCACAATCATGAAGCATTTTGACATATGTCGTTCCTGAGTTTTGTCCTGTCGCATCATAATGGAAATCCTAGGGCAAACGGGCTTTCGGTGCCA"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, pA[1:50])
v$viterbi_seq
template
## didnt seem to help


## we could potentially try a clustering approach for base-calling
## all the states are defined as cluster centers and we just see which points cluster to which states
##--> just imputing
## or we could do kmeans clustering where we assume there are 1024 states -- or do BIC to pick something nearby


## perhaps the HMM should be looking at 2 emissions: the mean AND the sd
## 1/16/14 -- just made 2emit versions of F, B, and V
########### TWO EMITS
source("~/searchPaths/bitbucket/squiggler/example_dtw/hmm_functions.R")
## made function to automatically generate transition probs
## put dimers, trimers, fourmers, and fivemers in environment
#Dimers:
dimer.t <- generate.kmer.transition.probs(dimers)
dimer.e.level <- generate.kmer.emission.probs(dimers)
dimer.e.sd <- generate.kmer.emission.probs(dimers, level=FALSE)
dimer.i <- generate.kmer.initial.probs(dimers)
answer <- generate.emissions.twoEmits(emissions1 = dimer.e.level, emissions2 = dimer.e.sd, transitions = dimer.t, initial = dimer.i, states = dimers, length=500) 

fl <- forward.long.twoEmits(emissions1 = dimer.e.level, emissions2 = dimer.e.sd, transitions = dimer.t, initial = dimer.i, states = dimers, emitted.data1 = answer$emitted.data1, emitted.data2 = answer$emitted.data2)
dl <- prob.data.long(Forward = fl$forward, scalefactors = fl$scales)
bl <- backward.long.twoEmits(emissions1 = dimer.e.level, emissions2 = dimer.e.sd, transitions = dimer.t, initial = dimer.i, states = dimers, emitted.data1 = answer$emitted.data1, emitted.data2 = answer$emitted.data2)
v <- viterbi.twoEmits(emissions1 = dimer.e.level, emissions2 = dimer.e.sd, transitions = dimer.t, initial = dimer.i, states = dimers, emitted.data1 = answer$emitted.data1, emitted.data2 = answer$emitted.data2)
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states = dimers)
compare.statepath(v$viterbi_path, answer$statepath)
compare.statepath(cl$centroid.path, answer$statepath)   
compare.statepath(cl$centroid.path, v$viterbi_path)
#compare seqs
compare.seq(v$viterbi_seq, answer$statepathseq)
compare.seq(cl$centroid.seq, answer$statepathseq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)
### CONCLUSION: TWO EMISSIONS DEF BETTER THAN ONE. MORE IS MORE. Centroid seems to perform better when they do not perform the same (which is often).
### ...did not seem to improve real data (did not do worse either)


## TWO EMITS -- try on real data
## SCIARA1 - TWO EMITS 54.7% ID with ONT -- HMM1=54.7%, 
fivemer.t <- generate.kmer.transition.probs(fivemers)
model_r7.3 <- read.model.file("template_model_r7.3.tsv")
fivemer.e.1 <- t(model_r7.3[,2:3])
fivemer.e.2 <- t(model_r7.3[,4:5])
data <- read.events.file("../data/events/sciara_ch111_file162_template.tsv")
emitted.data1 <- data$mean 
emitted.data2 <- data$stddev
sci_ch111_template="CTGGGTGTTTATACCAACGTGTTCTACAACGCATTATGCAGTATAATCCTCCATTTCTTAGTCGCCTTTATTACTAGAAATGTTTTCTAGTGAGGGATTCTTTCGGGCCCATGCCGAAATGAAACCAGCATTTCAAAGCGCTACCCCGCTCATGCTCGTATTTGTCCCGTGCATA"
v <- viterbi.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
fl <- forward.long.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
bl <- backward.long.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states = fivemers)
compare.seq(substr(sci_ch111_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(sci_ch111_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)

##SCIARA2 - TWO EMITS 72.2% with ONT -- HMM1=72.2%, 
data <- read.events.file("../data/events/sciara_ch224_file20_template.tsv")
emitted.data1 <- data$mean 
emitted.data2 <- data$stddev
sci_ch224_template="AGTATAGTCTCATCGGGTTAAGAAAATATTCTCGCAGGTTAGTTGTGTTTCTGTTAATACCATTTTACTACGACTTCGAGCCGTGAGGCCCGGGTGAGCTGGGCACAAGGATAGGAGCGACGTTCCTCTTCCACCTTAATCCTCCAGTGGTAGCCATTAGCCATGATGTAACAAT"
v <- viterbi.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
fl <- forward.long.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
bl <- backward.long.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states = fivemers)
compare.seq(substr(sci_ch224_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(sci_ch224_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)

##LAMBDA - TWO EMITS 57.1% ID with ONT -- HMM1=60.7%, 
model_r7 <- read.model.file("template_model_r7.tsv")
# fivemer.e <- t(model_r7[,2:3])
fivemer.e <- t(model_r7.3[,2:3]) ## r7.3 actually works better
data <- read.events.file("../data/events/lambda_ch230_file5_template.tsv")
emitted.data1 <- data$mean 
emitted.data2 <- data$stddev
lambda_ch230_template="TTTTCAATTTATGAAGGTCCGTCAGGCGTTCGTTCTTCTTCGTCATATCCTGGTTTCTTCATGAGCTCAGGGAAGGAGCACAATCATGAAGCATTTTGACATATGTCGTTCCTGAGTTTTGTCCTGTCGCATCATAATGGAAATCCTAGGGCAAACGGGCTTTCGGTGCCAATAT"
v <- viterbi.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
fl <- forward.long.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
bl <- backward.long.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states = fivemers)
compare.seq(substr(lambda_ch230_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(lambda_ch230_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)
v$moves
data$move[1:50]


#LAMBDA2 - TWO EMITS 48.1% ID with ONT -- HMM1=53.7%
data <- read.events.file("../data/events/lambda_ch304_file2_template.tsv")
emitted.data1 <- data$mean 
emitted.data2 <- data$stddev
lambda_ch304_template <- "GTGTTATGCTGATCCCATTATGGAAGTTTCATTAGAGACAAAGCCGTCGGGACTATCATTATCGATACTGAACTCCCAACGTGCGTGGAGCATCAACCGTTACGGCTTAATTATGAAAATACTACGTATTAATTATTATCCATCTATCTAACGTAAAATAGGTGCGACATCAAGT"
v <- viterbi.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
fl <- forward.long.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
bl <- backward.long.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, states = fivemers)
compare.seq(substr(lambda_ch304_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(lambda_ch304_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)



# messing around with these numbers....
# scale=0.896359165058
# shift = -1.83524398403
# drift = -0.000663759489
# scale_sd = 0.995414047068
# var_sd = 2.906893408
# var = 1.44664928588
# fivemer.e.1 <- t(model_r7.3[,2:3])
# fivemer.e.2 <- t(model_r7.3[,4:5])
# data <- read.table("../data/events/sciara_ch111_file162_template.tsv")
# emitted.data1 <- (data$V1/scale)-shift ## actually this is 2D data
# emitted.data2 <- data$V2/scale_sd
# sci_ch111_template="ACTCTCTGAGAACCGCACGGATCAGACTTGCGACGAGCCCGGCTCACCCAATTGCAGTAAACTCTGAGAAACAATGGCTTAGATCGCCGAAGTGAGGCCACGGCCGGGGTAGAACGAAACGGCTGAAAAGAGAAGAGAGATGAGAAAATTTCAAGCAGAGCACGGCGG"
# v <- viterbi.twoEmits(emissions1 = fivemer.e.1, emissions2 = fivemer.e.2, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data1[1:50], emitted.data2[1:50])
# v$viterbi_seq
# sci_ch111_template


### Need to allow it to transition to 5-mers that are not nearest neighbor, but 2 and maybe 3-5 away
### For the 2 door down neighbor approach alone, can just include all those kmers as possible transitions
###  e.g. instead of looking for only overlap k-1 mers, define TPs for overlapping k-2 mers
###     --> can still give uniform transition prob to all 
###     --> can then identify the moves when reading the states into a sequence
###     ---> if overlap(k-1mers) move = 1, add 1 to seq, elif overlap(k-2mers), move = 2, add 2 to seq.
### May also need to allow the possibility that a kmer emitted more than 2 events (move = 0)
### In future perhaps also modeling in duration as a 3rd emission would be useful.
### Even if as a post-processing step looking at homo-5mers and seeing if they should be expanded
### As these alt moves become allowable, it might also be useful to do the book-keeping of what
###  state(kmer)-call corresponds to what event number
### Interestingly -- in move=0 scenarios, this model will usually call homopolymers as transitions to themselves 
##      whereas it would say move=0 for any hetero-polymer.

## 1. make generate.kmer.transition.probs.withgaps
## 2. make something that can translate to sequence with gaps/inserts allowed... not all kmers will overlap by k-1
## 3. make function to generate seqs/state paths
source("~/searchPaths/bitbucket/squiggler/example_dtw/hmm_functions.R")

## trimers
trimer.t <- generate.kmer.transition.probs.withgaps(trimers)
trimer.e <- generate.kmer.emission.probs(trimers)
trimer.i <- generate.kmer.initial.probs(trimers)
answer <- generate.emissions(emissions = trimer.e, transitions = trimer.t, initial = trimer.i, states = trimers, length=50, withGaps = TRUE) 


fl <- forward.long(emissions = trimer.e, transitions = trimer.t, initial = trimer.i, states = trimers, answer$emitted.data)
dl <- prob.data.long(Forward = fl$forward, scalefactors = fl$scales)
bl <- backward.long(emissions = trimer.e, transitions = trimer.t, initial = trimer.i, states = trimers, answer$emitted.data)
v <- viterbi(emissions = trimer.e, transitions = trimer.t, initial = trimer.i, states = trimers, answer$emitted.data, getseqfxn = get.sequence.withgaps)
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, getseqfxn = get.sequence.withgaps, states = trimers)
compare.statepath(v$viterbi_path, answer$statepath)
compare.statepath(cl$centroid.path, answer$statepath)   
compare.statepath(cl$centroid.path, v$viterbi_path)

#compare seqs
compare.seq(v$viterbi_seq, answer$statepathseq)
compare.seq(cl$centroid.seq, answer$statepathseq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)
## STATES 80% ID ... SEQ 57-58% ID.... somehow putting seq together wrong...
### -- this is b/c need to algn needleman-wunsch style

## note to self -- scaling events might have to do with max(HairPin)/max(leader) OR max(HairPin)-max(leader) 
answer$moves
v$moves
cl$moves
answer$statepathseq
v$viterbi_seq
cl$centroid.seq


## Try GAPPED out on real data -- no bias correction
# run.gapped.real.data()
## SCIARA1 -- 43.1% ID with ONT; TWO EMITS=54.7%; HMM1=54.7%,;  
fivemer.t <- generate.kmer.transition.probs.withgaps(fivemers)
model_r7.3 <- read.model.file("template_model_r7.3.tsv")
fivemer.e <- t(model_r7.3[,2:3])
data <- read.table("../data/events/sciara_ch111_file162_template.tsv")
emitted.data <- data$V1 
sci_ch111_template="CTGGGTGTTTATACCAACGTGTTCTACAACGCATTATGCAGTATAATCCTCCATTTCTTAGTCGCCTTTATTACTAGAAATGTTTTCTAGTGAGGGATTCTTTCGGGCCCATGCCGAAATGAAACCAGCATTTCAAAGCGCTACCCCGCTCATGCTCGTATTTGTCCCGTGCATA"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50], getseqfxn = get.sequence.withgaps)
fl <- forward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
bl <- backward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, getseqfxn = get.sequence.withgaps, states = fivemers)
compare.seq(substr(sci_ch111_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(sci_ch111_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)
# Try on input events -- where I can see that they do not equal the template events reported (prob due to 0 moves that cannot be connected to a new kmer)
inputdata <- read.events.file(events.file = "../data/events/sciara_ch111_file162_input.tsv",input.events = TRUE)
emitted.input <- inputdata$mean[51:(51+length(emitted.data)-1)]
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.input[1:50], getseqfxn = get.sequence.withgaps)
fl <- forward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.input[1:50])
bl <- backward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.input[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, getseqfxn = get.sequence.withgaps, states = fivemers)
compare.seq(substr(sci_ch224_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(sci_ch224_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)





##SCIARA2 -- 73.58% ID with ONT; TWO.EMITS=72.2%; HMM1=72.2%,; 
fivemer.t <- generate.kmer.transition.probs.withgaps(fivemers)
data <- read.table("../data/events/sciara_ch224_file20_template.tsv")
emitted.data <- data$V1 
sci_ch224_template="AGTATAGTCTCATCGGGTTAAGAAAATATTCTCGCAGGTTAGTTGTGTTTCTGTTAATACCATTTTACTACGACTTCGAGCCGTGAGGCCCGGGTGAGCTGGGCACAAGGATAGGAGCGACGTTCCTCTTCCACCTTAATCCTCCAGTGGTAGCCATTAGCCATGATGTAACAAT"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50], getseqfxn = get.sequence.withgaps)
fl <- forward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
bl <- backward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, getseqfxn = get.sequence.withgaps, states = fivemers)
compare.seq(substr(sci_ch224_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(sci_ch224_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)
# Try on input events -- where I can see that they do not equal the template events reported (prob due to 0 moves that cannot be connected to a new kmer)
inputdata <- read.events.file(events.file = "../data/events/sciara_ch224_file20_input.tsv",input.events = TRUE)
emitted.input <- inputdata$mean[51:(51+length(emitted.data)-1)]
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.input[1:50], getseqfxn = get.sequence.withgaps)
fl <- forward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.input[1:50])
bl <- backward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.input[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, getseqfxn = get.sequence.withgaps, states = fivemers)
compare.seq(substr(sci_ch224_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(sci_ch224_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)

##LAMBDA -- 78.3% .. or 72.13% ID with ONT; TWO.EMITS=57.1% ; HMM1=60.7%, ; 
model_r7 <- read.model.file("template_model_r7.tsv")
fivemer.e <- t(model_r7[,2:3])
fivemer.e <- t(model_r7.3[,2:3]) 
data <- read.table("../data/events/lambda_ch230_file5_template.tsv")
emitted.data <- data$V1 
lambda_ch230_template="TTTTCAATTTATGAAGGTCCGTCAGGCGTTCGTTCTTCTTCGTCATATCCTGGTTTCTTCATGAGCTCAGGGAAGGAGCACAATCATGAAGCATTTTGACATATGTCGTTCCTGAGTTTTGTCCTGTCGCATCATAATGGAAATCCTAGGGCAAACGGGCTTTCGGTGCCAATAT"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50], getseqfxn = get.sequence.withgaps)
compare.seq(substr(lambda_ch230_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
v$moves
data$V7[1:50]
#150bp
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:150], getseqfxn = get.sequence.withgaps)
compare.seq(substr(lambda_ch230_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
v$moves
data$V7[1:length(v$moves)]

#LAMBDA2 -- 55% ID with ONT; TWO.EMITS=48.1% ; HMM1=53.7% ;
data <- read.events.file("../data/events/lambda_ch304_file2_template.tsv")
emitted.data <- data$mean
lambda_ch304_template <- "GTGTTATGCTGATCCCATTATGGAAGTTTCATTAGAGACAAAGCCGTCGGGACTATCATTATCGATACTGAACTCCCAACGTGCGTGGAGCATCAACCGTTACGGCTTAATTATGAAAATACTACGTATTAATTATTATCCATCTATCTAACGTAAAATAGGTGCGACATCAAGT"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50], getseqfxn = get.sequence.withgaps)
compare.seq(substr(lambda_ch304_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
v$moves
data$move[1:50]
#150bp
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:150], getseqfxn = get.sequence.withgaps)
compare.seq(substr(lambda_ch304_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)

#LAMBDA3 -- % ID with ONT; TWO.EMITS=__ ; HMM1=__ ;
fivemer.e <- t(model_r7[,2:3])
fivemer.e <- t(model_r7.3[,2:3])
data <- read.events.file("../data/events/lambda_ch404_file4_template.tsv")
emitted.data <- data$mean
lambda_ch404_template <- "GCCTCAACACAGATGTCCCTATTTAGCAAGGTTTTTCCGTTAAGGCCGTTTCCGTCTTCTTCGTCACTGTAATGTTTTATTTAAATATAGACGGCGGAAAGGAAGGCACAGATCGGGCAAGCATTTTGGCCTCTGTCGTTCATTTCTCTGTTTTGTCCGTGGAATGAACAATGGAAGTTAGCTCAAGGCAGCTGAACCCCTTTCGGTGCGTATAGCCGTACCGGTCCAGAGGATAACTGGAGGAGCCGGAATCCCGTTCTGCGGCAATGGGGGTAATGAGGTGCTTTATGATGCCGCCTTCAACACTAAAATGATATGCCGAAAGGGATGCTGGCATTGAAGAACGAAGG"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50], getseqfxn = get.sequence.withgaps)
compare.seq(substr(lambda_ch404_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
v$moves
data$move[1:50]
#150bp
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:150], getseqfxn = get.sequence.withgaps)
compare.seq(substr(lambda_ch404_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
## whole thing -- r7 took 84 minuts, r7.3 took 70 minutes
date()
fivemer.e <- t(model_r7[,2:3])
vr7 <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data, getseqfxn = get.sequence.withgaps)
date()
fivemer.e <- t(model_r7.3[,2:3])
vr73 <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data, getseqfxn = get.sequence.withgaps)
date()
## compare to lambda genome
lambdaGenome <- readDNAStringSet(filepath = "/Users/johnurban/data/otherGenomes/Lambda/lambdaGenomeSequence.fa", format = "fasta", use.names = FALSE)
ont <- readDNAStringSet(filepath = "/Users/johnurban/searchPaths/bitbucket/squiggler/data/fasta/lambda_ch404_file4.fa", format = "fasta", use.names = FALSE)
compare.seq(substr(lambdaGenome,1,2000), substr(vr7$viterbi_seq,1,2000)) # 1:2000, 56.8% ID
compare.seq(substr(lambdaGenome,1,2000), substr(vr73$viterbi_seq,1,2000)) # 1:2000, 45.1% ID
compare.seq(substr(lambdaGenome,1,2000), substr(ont,1,2000)) # 1:2000, 66.6% ID
compare.seq(substr(vr7$viterbi_seq,1,2000), substr(ont,1,2000)) # 1:2000, 65.4% ID
compare.seq(substr(vr73$viterbi_seq,1,2000), substr(ont,1,2000)) # 1:2000, 46.3% ID


##L4 -- 20150225

fivemer.t <- generate.kmer.transition.probs(fivemers)
fivemer.e <- t(model_r7[,2:3])
fivemer.e <- t(model_r7.3[,2:3])
fivemer.i <- generate.kmer.initial.probs(fivemers)
data <- read.table("../data/events/lambda_ch204_file2_notbasecalled_template_events.tsv")
emitted.data <- data$V1 
templateONT <- "AGTAAGCTCCAATGTTGGCCCTTAGCGCAGCATTCAGTACTTGCATGATTCATAGGTCTCTATGCCGTGAGATACTTTACCCGTGATAGTCTCAGTGGCGAGAAGTAATATCGATCTTTTCGTCATAAATTAGGGGATAATCGCACAAACGAACGTAGTTACTAATTCTTCCCGGTCGTAGTCATAGAAAAATTCGACGTACATTCGGTGATTATGTCGCTAAATGGGTTGTCGTGGAATAGCCCCAATATATACTCGTAATTCTAGTACTGCTCAGAAGCGGGGTAGAATAGCTGCGTGGGGATGGGATAGGCGCGTTGCAGCGTGGTCGTACGGGCAGTCCAGTTTCAAAGCTGGGGTAAATGTCCCATGAACGCAAGTGGGGAGATCACGAAGAGAGAGATAGCCCCCATGAACAGACGGCACCCGGCATTGGTTTCTTTCGGATCCGTAAGTTCTTCACGTTGATTGTGCCAAAGTCTTATCCGGTTAGTATGCACCGCTGATGAGTGGACTTGGCAGTCCTTTGATGTCTTTTGTCTAGCCGTGCTAATAACGTTGAATGACCCAGGGCATTCGTCAGGTTCTTCATGTATCTGGAGGTCTTATGGCATTCTCCTGACCCCCGACACGCTCGAGAATAACCAGAGAAATGGAGCGTTTTCGAGAAAGCAGCTACAATCTCGAAGTGATCTGCAGAATTGTATTCTTAATTCTTGAGTTGTTAATCATCGCACTATGTTTTGTATCCAGAGAGACGAAAGAACATCGTTGAGCCCTAGTCAATCTACATAGTCTTTGTGCACGTGGGCATTAGAACGTCACCAGAGGATTGCGAATGGGAGGGATTGTCTTCATGCTATAGTACCCACTCTTCGAAATTTATCTTATGCCGCAGATACATATCACCACAACTATTCATCTACATAAGCGCTTCAACGTACTTCTGCATGGAAGGCGCTAAGATTATCTGCACATCAAATCTATAGCTAGAGGAATAATCTTTTTAACGGGATACCTGGATCTTGATACCATTTATAATATATAACTCATTGTGCGAAGAGCCTCTTATGTCAGTTCATAGAATGCTGCTGGCGGCAGAAATGTTATGTGAAGCGGCTCAGATAACACTCCGGGCGATGTTGTGAACCTCATTCGACATGGTAGGACGGAGAGATATGATAACACTGAAACAATCGGCATAAGCTCGGACGCGGTTCCAGGATTTCACCAGCGCGTCTGTTTACAATGAATCTACACGCCCTACAGCTCGCTTTGATGCGAATCTGCGCGCGCGCGCTTCCAGCACCGACCACCCTATGCAATACCCTCGCTGCATCCTCTATCGCATTATATACAAACCCCAGAATCAAGCGCGGCTTATTGGAAGTAATGGTCTTCACATAGTTACTAAAATAGGGTAGGTCTTGTTCGATTGGTCACTATCCGCAGGGGAAGACCTTCCATTACGGCAGTTGTGGAGAGTCGTACTACCGCGAGGAATCTTCCCGAAGGATGATGAAATAAGAAACTATACGGTTTTTCTTTTCGTCGTGCTGCCAGAGACAACTTCGTCCCACTTATTTCCGTTGGCAATACTACGGTAACTCGAGATACGAAGATAACACCCACCCAGCGTTCTGTTTAG"
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50], getseqfxn = get.sequence.withgaps)
fl <- forward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
bl <- backward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.data[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, getseqfxn = get.sequence.withgaps, states = fivemers)
compare.seq(substr(templateONT,1,50), v$viterbi_seq)
compare.seq(substr(templateONT,1,50), cl$centroid.seq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)
# Try on input events -- where I can see that they do not equal the template events reported (prob due to 0 moves that cannot be connected to a new kmer)
inputdata <- read.events.file(events.file = "../data/events/sciara_ch224_file20_input.tsv",input.events = TRUE)
emitted.input <- inputdata$mean[51:(51+length(emitted.data)-1)]
v <- viterbi(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.input[1:50], getseqfxn = get.sequence.withgaps)
fl <- forward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.input[1:50])
bl <- backward.long(emissions = fivemer.e, transitions = fivemer.t, initial = fivemer.i, states = fivemers, emitted.input[1:50])
cl <- centroid.long(Forward = fl$forward, F_scales = fl$scales, Backward = bl$backward, B_scales = bl$scales, getseqfxn = get.sequence.withgaps, states = fivemers)
compare.seq(substr(sci_ch224_template,1,nchar(v$viterbi_seq)), v$viterbi_seq)
compare.seq(substr(sci_ch224_template,1,nchar(cl$centroid.seq)), cl$centroid.seq)  
compare.seq(cl$centroid.seq, v$viterbi_seq)
###########################################################
## alignments from bioconductor biostrings
# s1 <- "GGGAGCAAGCTACCGCGCCAATAGGGAATAAGAACAACTGCCAAGTGAGCATCA"
# s2 <- "GGGAGCAAGCTACTCGCGCCAATAGGGAATAAGAACAGTCGCCAAGTGAGCATCA"
# ## minimize mismatches and gaps
# sigma <- nucleotideSubstitutionMatrix(match = 0, mismatch = -1, baseOnly = TRUE)
# pairwiseAlignment(s1, s2, substitutionMatrix = sigma, gapOpening = -1, gapExtension = -1, scoreOnly = FALSE)
# # maximize matches
# sigma <- nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = TRUE)
# pairwiseAlignment(s1, s2, substitutionMatrix = sigma, gapOpening = 0, gapExtension = 0, scoreOnly = FALSE)
###########################################################


######### GAPS AND TWO EMITS ##############################














#############OTHER NOTES########
## can try as a pair HMM
## where there are 3 matrices:
##  Match -- where emission matches state
##  Gap in emissions -- where there is only a transition prob to the gap, no emission (i.e. skip a 5mer, move 2)
##  Gap in data -- where there is an orphaned emission with only a transition probability, no new data point (i.e. move = 0)
## An advantage of this is...
## This will directly align emissions to kmers.
###   e1 e2 -- e3 e4 e5 e6
###   k1 k2 k3 k4 k5 -- k6
### I guess the problem still is one of defining allowlable transitions

## remmeber to change functions s.t. the transition probs are logProbs in the matrix so this operation does not have to be done every time
## also it might be faster to just store the locations of possible transitions given suffixes as I have done for other things
##  -- e.g. for the viterbi, you know all transitions with prob=0 will NOT be max, so just look at those that might be 
##        instead of looping over all 1024 states every time
