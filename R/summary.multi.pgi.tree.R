summary.multi.pgi.tree <-
function(object) {
	cat(length(object)," PGi Data objects\n")
	status.t<-sapply(object,function(x){return(x$type)})
	cat("Status of the objects: ",paste(status.t,collapse=","),"\n")
	if(any(status.t == "filled") || any(status.t == "consensus"))	cat("Tree lengths recorded in the PGi Objects: ",paste(sapply(object,function(x){return(x$tree.length)}),collapse=","),"\n")
}

