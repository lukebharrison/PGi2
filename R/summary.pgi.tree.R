summary.pgi.tree <-
function(object,print.all.seqs=FALSE) {
	cat("**PGi data object**\n")
	if(object$type=="empty") {
		cat("Contains a phylogenetic tree topology and sequence data\n")
	} else if(object$type == "filled") {
		cat("Contains a filled phylogenetic tree\n")
	} else if(object$type == "consensus") {
		cat("Contains a PGi Consensus tree\n")
	} else if(object$type == "partial") {
		cat("Contains a partially filled phylogenetic tree\n")
	} else if(object$type == "superconsensus") {
		cat("Contains a PGi superconsensus tree of multiple pseudoconseus trees\n")
	}
	cat("Number of taxa: ",object$ntaxa,"\n",sep="")
	cat("Number of nodes: ",object$nnodes,"\n",sep="")
	cat("Sequence Length: ",object$seq.len,"\n",sep="")
	cat("\n*Topology*\n")
	t.obj<-object
	t.obj$type<-"empty"
	t.phylo<-pgi.to.phylo(t.obj)
	cat(write.tree(t.phylo))
	cat("\n\n")
	cat("*Data* (First two sequences only, override with print.all.seqs=TRUE)\n")
	if(print.all.seqs) max<-object$ntaxa else max<-min(2,object$ntaxa)
	align<-max(sapply(t.phylo$tip.label,nchar))+5
	for(i in 1:max) cat(t.phylo$tip.label[i],":",rep(" ",times=max(0,align-nchar(t.phylo$tip.label[i]))),t.phylo$tip.seqs[i],"\n",sep="")
	cat("\n")

	if(object$type=="partial") {
		cat("*Status of Partial PGi Execution\n")
		cat(.pgi.getn.completed.nodes(pgi.tree),"nodes out of",pgi.tree$nnodes,"filled\n")
	}
	if(any(object$type==c("filled","consensus"))) {
		cat("*PGi Execution Details*\n")
		cat("Tree length: \t\t",object$tree.length," sequence heterochronies\n\n",sep="")
		cat("PGi was executed with the following parameters:\n")
	} else
	if((object$type == "partial" || object$type == "empty") && (!is.null(object$inf.params))) {
		cat("PGi execution parameters:\n")
	}
	if(!is.null(object$inf.params)) {
		align<-max(sapply(names(object$inf.params),nchar))+5
		for(i in 1:length(object$inf.params)) {
			if(names(object$inf.params)[i] != "edit.cost.mat" || names(object$inf.params)[i] != "all.dev.seq.mat") {
				cat(names(object$inf.params)[i],":",rep(" ",times=max(0,align-nchar(names(object$inf.params)[i]))),object$inf.params[[i]],"\n",sep="")
			}
		}
		cat("\n")
	} 
	if(object$type=="consensus") {
		cat("*PGi Pseudoconsensus Details*\n")
		cat("The PGi Peudoconsensus was calculated with the following parameters:\n")
	} else if(any(object$type == c("empty","filled")) && (!is.null(object$inf.params))) {
		cat("PGi Consensus parameters:\n")
	}
	if(!is.null(object$con.params)) {
		align<-max(sapply(names(object$con.params),nchar))+5
		for(i in 1:length(object$con.params)) cat(names(object$con.params)[i],":",rep(" ",times=max(0,align-nchar(names(object$con.params)[i]))),object$con.params[[i]],"\n",sep="")
		cat("\n")
	}

	## Superconesnus
	if(object$type=="superconsensus") {
		cat("*PGi Superconsensus Details*\n")
		cat("The PGi superconsensus was calculated with the following parameters:\n")
		cat("Number of Pseudoconsensus Trees Inputed: ",length(object$supercon.TreeLengths),"\n",sep="")
		cat("Tree Lengths of Inputed Pseudoconsensus Trees: ",paste(object$supercon.TreeLengths,collapse=", "),"\n",sep="")
		if(object$supercon.tol > 0) {
			cat("Maximum Allowable Deviation from Minimum Tree Length: ",object$supercon.tol,"%, or ",(min(object$supercon.TreeLengths)*object$supercon.tol)," sequence heterochronies\n",sep="")
			cutoff<-min(object$supercon.TreeLengths)+(min(object$supercon.TreeLengths)*object$supercon.tol)
			nused.trees<-length(which(object$supercon.TreeLengths <= cutoff))
			cat("Superconsensus was calculated using ",nused.trees," Pseudoconesus trees below the cutoff\n",sep="")
		} else {
			cat("Superconsesus was calculated by explicitely using ALL input pseudoconsensus trees, regardless of tree length\n")
		}
	}
	cat("\n")

}


