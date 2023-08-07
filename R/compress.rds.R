compress.rds <-
function(rds) {
 offset<-0
 new_rds<-rds
 sel_rds<-rds
 for(i in 1:length(unique(rds))) {
  min<-min(sel_rds)
  sel_rds[sel_rds == min]<-999
  if(min == -2) {
   new_rds[new_rds == min] <- -2
   offset<-offset+1
  } 
  else {
   new_rds[new_rds == min]<-(i-offset)
  }
 }
 return(new_rds)
}

