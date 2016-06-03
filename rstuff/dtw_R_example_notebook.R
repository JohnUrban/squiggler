## set pwd to ~/path/to/squiggler/example_dtw/ 
setwd("~/searchPaths/bitbucket/squiggler/example_dtw/") 

loadEvents <- function(eventsfile){
  data <- read.table(eventsfile)
  colnames(data) <- c("mean", "stddev", "start", "length")
  return(data)
}

sc1.in <- loadEvents("../data/events/sciara_ch111_file162_input.tsv")
sc2.in <- loadEvents("../data/events/sciara_ch224_file20_input.tsv")
l1.in <- loadEvents("../data/events/lambda_ch230_file5_input.tsv")
l2.in <- loadEvents("../data/events/lambda_ch304_file2_input.tsv")

## params
ylim <- range(sc1.in$mean, sc2.in$mean, l1.in$mean, l2.in$mean)

## example plot 1
plot(1:50, sc1.in$mean[1:50], type="l", col="dark blue", ylim=ylim, ylab="event height", xlab="event number")
lines(1:50, sc2.in$mean[1:50], col="dark cyan")
lines(1:50, l1.in$mean[1:50], col="dark red")
lines(1:50, l2.in$mean[1:50], col="dark orange")

## example plot 2
plot(sc1.in$start[1:50]-sc1.in$start[1], sc1.in$mean[1:50], type="l", col="dark blue", ylim=ylim, ylab="event height", xlab="relative time")
lines(sc2.in$start[1:50]-sc2.in$start[1], sc2.in$mean[1:50], col="dark cyan")
lines(l1.in$start[1:50]-l1.in$start[1], l1.in$mean[1:50], col="dark red")
lines(l2.in$start[1:50]-l2.in$start[1], l2.in$mean[1:50], col="dark orange")


## example dtw
library(dtw)
#alignment<-dtw(query,template,keep=TRUE)
alignment<-dtw(sc1.in$mean[1:50],sc2.in$mean[1:50],keep=TRUE)
plot(alignment,type="threeway")
plot(alignment,type="twoway")
alignment2 <- dtw(sc1.in$mean[1:50],sc2.in$mean[1:50],keep=TRUE, step=rabinerJuangStepPattern(6,"c"))
plot(alignment2, type="twoway",offset=-2);

## example dtw -- lambda
alignment<-dtw(l1.in$mean[1:50],l2.in$mean[1:50],keep=TRUE)
plot(alignment,type="threeway")
plot(alignment,type="twoway")
alignment2 <- dtw(sc1.in$mean[1:50],sc2.in$mean[1:50],keep=TRUE, step=rabinerJuangStepPattern(6,"c"))
plot(alignment2, type="twoway",offset=-2);




## See the recursion relation, as a figure and text
plot(rabinerJuangStepPattern(6,"c"))
rabinerJuangStepPattern(6,"c")
