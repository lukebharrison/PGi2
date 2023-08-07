seq.to.rds <-
function(dev_seq,seq_len) {
 offset<-0
 rds_dev_seq<-rep(0,times=seq_len)
 for(i in 1:length(dev_seq)) {
  for(j in 1:length(dev_seq[[i]])) {
   if(dev_seq[[i]][j] == -1) {
    if(length(dev_seq[[i]]) < 2) {
     offset <- offset + 1
    }
    next
   }
   rds_dev_seq[dev_seq[[i]][j]]<-(i-offset)
  }
 }
 rds_dev_seq<-rds_dev_seq[1:seq_len]
 rds_dev_seq[rds_dev_seq == 0] <- -2
 return(rds_dev_seq)
}

