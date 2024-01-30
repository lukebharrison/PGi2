pgi.inference <-
function(pgi.tree, verbosity=1, inf.params=NULL) {
	if(is.null(inf.params)) {
		if(is.null(pgi.tree$inf.params)) stop("No inference parameters specified in the tree or explicitely. Stopping")
	} else {
		cat("Overiding inference parametesr in the tree with provided parameters\n")
		pgi.tree$inf.params<-inf.params
	}
	if(pgi.tree$inf.params$heuristic == "sw04") {
		if(is.null(pgi.tree$inf.params$edit.cost.matrix)) {
			editcostmat<-create.edit.cost.matrix(pgi.tree,verbosity=verbosity)
			pgi.tree$inf.params$all.dev.seq.mat<-editcostmat[[1]]
			pgi.tree$inf.params$edit.cost.matrix<-editcostmat[[2]]
		}
	}
	if(pgi.tree$inf.params$heuristic == "exhaustive") {
		if(is.null(pgi.tree$inf.params$all.dev.seq.mat)) {
			pgi.tree$inf.params$all.dev.seq.mat<-unique(rbind(create.all.dev.seq.mat(pgi.tree$seq.len,pgi.tree$inf.params$simultaneity,verbosity=verbosity),extract.seqs.from.tree(pgi.tree$tree,pgi.tree$seq.len,c())))
		}
	}
	if(verbosity > 0) { 
		if(any(pgi.tree$inf.params$heuristic == c("pgi","oldgenetic"))) {
			cat("Running inference with the ",pgi.tree$inf.params$heuristic," heuristic. Replicates/cycle=",pgi.tree$inf.params$replicates," Cycles=",pgi.tree$inf.params$cycles," and a retained matrix size of ",pgi.tree$inf.params$ret.anc.seq,"\n",sep="") 
		} else {
			cat("Running inference with the ",pgi.tree$inf.params$heuristic," algorithm\n",sep="")
		}
	}
	if(verbosity == 1) cat("Will display overall (per run) progress bar\n")
	if(verbosity == 2) cat("Will display PER node progress bars. Note: there are ",pgi.tree$nnodes," nodes in total\n",sep="")
	if(verbosity == 1) {
		pb <- txtProgressBar(min=0, max = pgi.tree$nnodes, style = 3)
		setTxtProgressBar(pb,.pgi.getn.completed.nodes(pgi.tree))
	}	else pb<-NULL
	if(!pgi.tree$type == "partial") {
		pgi.tree$partial.tree<-pgi.tree$tree
		pgi.tree$type<-"partial"
	}
	pgi.tree$filled.tree<-.downsweep(pgi.tree,pgi.tree$partial.tree,progress=pb,verbosity=verbosity)[[2]]
	pgi.tree$type<-"filled"
	pgi.tree$partial.tree<-NULL
	if(pgi.tree$inf.params$resume) unlink(pgi.tree$inf.params$resume.temp.file)
	if(inf.params$use.simplecon.treelength) {
 	 t.simple.con<-pgi.pseudo.consensus(pgi.tree,con.params=list(con.type="simple",edit.cost.func=pgi.tree$inf.params$edit.cost.func),verbosity=verbosity)
	 pgi.tree$tree.length<-pgi.con.tree.length(t.simple.con$consensus)
	} else {
	 pgi.tree$tree.length<-min(pgi.tree$filled[[1]][,pgi.tree$seq.len+1])
	}
	if(verbosity == 1) close(pb)	
	if(verbosity > 2) summary(pgi.tree)
	return(pgi.tree) 
}

