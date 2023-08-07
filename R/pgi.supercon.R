pgi.supercon <-
function(con_trees,verbosity=1,tol=0) {
	if(class(con_trees) != "multi.pgi.tree") stop("Please provide a multi.pgi.tree object")
	if(!(all(sapply(con_trees,function(x){return(x$type)}) == "consensus"))) stop("All PGi data objects in the multi.pgi.tree object *MUST* have a calculated pseudoconsensus")
	tls<-	sapply(con_trees,function(x) { return(x$tree.length) })
	if(tol >= 0) {
		cutoff<-min(tls)+(min(tls)*tol)
		tree.idx<-which(tls <= cutoff)		
		if(verbosity > 0) cat("Calculating superconsensus using those trees with a tree length <= ", cutoff,"\n")
	} else {
		tree.idx<-seq(1,length(tls))
		if(verbosity > 0) cat("Calculating superconsensus using ALL trees\n")
	}

	if(verbosity > 1) cat("Selected",length(tree.idx),"trees. Pseudoconsensus tree numbers:",paste(tree.idx,collapse=","),"\n")
	supercontree<-con_trees[[tree.idx[1]]]
	cons.to.send<-c()
	supercontree$tree.length<-0
	for(i in tree.idx) {
		cons.to.send<-c(cons.to.send,list(con_trees[[i]]$consensus))
		supercontree$tree.length<-supercontree$tree.length+con_trees[[i]]$tree.length
	}
	supercontree$consensus<-.pgi.supercon.fill(supercontree$consensus,cons.to.send,length(tree.idx),supercontree$seq.len,verbosity=verbosity)
	supercontree$tree.length<-supercontree$tree.length/length(tree.idx)
	supercontree$filled.tree<-NULL
	supercontree$type="superconsensus"
	supercontree$con.params<-NULL
	supercontree$inf.params<-NULL
	supercontree$supercon.tol<-tol
	supercontree$supercon.TreeLengths<-tls
	if(verbosity > 1) cat("Done.\n")
	return(supercontree)
}

