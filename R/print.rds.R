print.rds <-
function(rds) {
 seq<-rds.to.seq(rds)
 if(length(seq[[1]])> 1) {
  nice_seq<-paste("[",gsub(" ","",toString(seq[[1]],trim=TRUE)),"]",sep="")
 } else {
  nice_seq<-paste(gsub(" ","",toString(seq[[1]],trim=TRUE)),sep="")
 }
	if(length(seq) > 1) {
	 for(i in 2:length(seq)) {
  	if(length(seq[[i]])> 1) {
  	 nice_seq<-paste(nice_seq,",[",gsub(" ","",toString(seq[[i]],trim=TRUE)),"]",sep="")
  	} else {
  	 nice_seq<-paste(nice_seq,",",gsub(" ","",toString(seq[[i]],trim=TRUE)),sep="")
  	}
	 }	
	}
 return(nice_seq)
}

