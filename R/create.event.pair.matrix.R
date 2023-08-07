create.event.pair.matrix <-
function(dev_seq) {
 event_pair <- array(0, dim=c(length(dev_seq)-1,length(dev_seq)-1))
 for (i in 2:length(dev_seq)-1) {
  for (j in (i+1):length(dev_seq)) {
   if(dev_seq[j] > dev_seq[i]) event_pair[j-1,i] <- 2
   if(dev_seq[j] < dev_seq[i]) event_pair[j-1,i] <- 0
   if(dev_seq[j] == dev_seq[i]) event_pair[j-1,i] <- 1
   if(dev_seq[j] == -2 || dev_seq[i] == -2) event_pair[j-1,i] <- "?"
  }
 }
 return(event_pair)
}

