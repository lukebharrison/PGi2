rds.to.seq <-
function(rds_dev_seq) {
 ln<-unique(as.vector(rds_dev_seq))[unique(as.vector(rds_dev_seq)) > 0]
 dev_seq<-as.list(rep(c(0),times=length(ln)))
 for(i in 1:length(rds_dev_seq)) {
  if(rds_dev_seq[i] == -2 || rds_dev_seq[i] == -1) {
   next
  }
  if(dev_seq[[rds_dev_seq[i]]][1] != 0) {
   dev_seq[[rds_dev_seq[i]]] <- c(dev_seq[[rds_dev_seq[i]]],i)
  } else {
   dev_seq[[rds_dev_seq[i]]] <- i
  }
 }
 return(dev_seq)
}

