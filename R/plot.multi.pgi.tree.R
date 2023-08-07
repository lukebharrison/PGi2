plot.multi.pgi.tree <-
function(multi.pgi.tree, format=NULL, ...) {
	## Function to plot trees
	## Format = rows / column to make a multiple plot or sequential in the case of c(1,1)
	num.trees<-length(multi.pgi.tree)
	if(is.null(format)) {
	 format.y<-floor(sqrt(num.trees))
	 format.x<-ceiling(num.trees/format.y)
	} else if(length(format) == 2) {
	 # presume a valid format command if correct length
	 format.y<-format[2]
	 format.x<-format[1]
	} else {
	 cat("Error, invalid format\n")
	 return()
	}
	cat("Creating multi.pgi.tree plot using a plotting matrix of x = ",format.x," and y = ",format.y,"\n",sep="")
	par(mfrow=c(format.x,format.y),ask=TRUE)
	for(i in 1:num.trees) {
		plot.pgi.tree(multi.pgi.tree[[i]],...)
	}
}

