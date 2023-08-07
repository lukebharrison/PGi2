pgi <-
function(pgi.tree, interactive=FALSE, nruns=NULL, verbosity = 0, replicates = 100, cycles = 100, ret.anc.seq = 100, simul = TRUE, heuristic ="pgi", edit.cost.func="parsimov", con.type="exhaustive", semi.ex.con.max.n=0, supercon=FALSE,supercon.tol=0,overwrite=FALSE, inf.params=NULL, con.params=NULL) {
	if(missing(pgi.tree) || (class(pgi.tree) != "pgi.tree")) {
		if(interactive) {
			cat("No PGi data object provided, attemping to load from a user-specified file\n")
			pgi.tree<-pgi.read.nexus(interactive=TRUE)
		} else {
			stop("PGi data object missing or not in pgi.tree format (use pgi.read.nexus() or pgi.assemble.tree()) and interactive loading is disabled")
		}
	}
	if(!pgi.tree$type == "empty") {
		if(!(overwrite || verbosity < 0)) stop("WARNING: This PGi data object has already been computed, continuing will overwrite exisitng runs/consensus. Stopping. To overwrite, add 'overwrite=TRUE'\n")
	}
	if(interactive) pgi.tree<-pgi.interactive.params(pgi.tree)
	if(!is.null(pgi.tree$verbosity)) verbosity <- pgi.tree$verbosity
	if(is.null(inf.params)) {
		if(is.null(pgi.tree$inf.params)) {
			if(verbosity >= 0) cat("WARNING: No inference parameters provided, using defaults. Or: inference parameters have been specifically specified to the PGi function indvidually\n")
			pgi.tree$inf.params<-list(replicates=replicates,cycles=cycles,ret.anc.seq=ret.anc.seq,simultaneity=simul,heuristic=heuristic,edit.cost.func=edit.cost.func,resume=FALSE,resume.temp.file="temp.execution.rda")
		}
	}	else {
		if(verbosity >= 0) cat("Overriding inf parameters in tree with the provided list\n")
		pgi.tree$inf.params<-inf.params
	}
	if(is.null(con.params)) {
		if(is.null(pgi.tree$con.params)) {
			if(verbosity >= 0) cat("WARNING: No consensus parameters provided, using defaults. Or: consensus parameters have been specifically specified to the PGi function indvidually\n")
			pgi.tree$con.params<-list(con.type=con.type,semi.ex.con.max.n=semi.ex.con.max.n,edit.cost.func=edit.cost.func)
		} 
	} else {
		if(verbosity >= 0) cat("Overriding con parameters in tree with the provided list\n")
		pgi.tree$con.params<-con.params
	}
	if(is.null(pgi.tree$inf.params$resume)) pgi.tree$inf.params$resume<-FALSE
	if(pgi.tree$inf.params$resume) {
		cat("Will save temporary file,",pgi.tree$inf.params$resume.temp.file,"with resuming informating if the computation is aborted\n")
		if(file.exists(pgi.tree$inf.params$resume.temp.file)) {
			cat("Found a matching temporary information file; will resume the computation if there is progress\n")
			# asked to resume and there is a resume file
			if(exists("partial.pgi.tree")) rm(partial.pgi.tree)
			load(pgi.tree$inf.params$resume.temp.file)
			if(partial.pgi.tree$type == "filled") {
				# make sure we account for a possible change in consensus parameters
				partial.pgi.tree$con.params<-pgi.tree$con.params
				pgi.tree<-partial.pgi.tree
			} else if(.pgi.getn.completed.nodes(partial.pgi.tree) > 0) {
				pgi.tree$partial.tree<-partial.pgi.tree$partial.tree
				pgi.tree$type<-"partial"
				rm("partial.pgi.tree")
			}
		}
	}
	if(is.null(nruns)) {
		if(is.null(pgi.tree$nruns)) {
			if(verbosity > 0) cat("WARNING: Nruns not defined explicitely or in the input PGi object, using c(1,1) as a default\n")
			nruns<-c(1,1)	
		} else {
			nruns<-pgi.tree$nruns	
		}
	} else {
		if(verbosity > 0) cat("Overiding nruns with the explicitely defined value\n")
	}
	total.runs<-nruns[1]*nruns[2]
	if((total.runs > 1) && (pgi.tree$type=="partial")) {
		stop("Error: can't resume a parallel/serial run; only a single run (for now)\n")
	}
	if(verbosity > 0) { 
		cat("Executing PGi, with ",total.runs," independent computation[s], total\n",sep="")
	}
	if(nruns[2] > 1) {
		parallel<-TRUE 
		require(foreach)
		require(doMC)
		registerDoMC(cores=nruns[2])	
		#init<-(foreach(j=1:nruns[2]) %dopar% (1+1))
		#if(verbosity > 2) {
		#	print("Init matrix for foreach")
		#	print(init)
		#}
		if(verbosity > 0) cat("Running PGi in parallel mode: running ",nruns[2]," simultaenous computations a total of ",nruns[1]," time[s].\nWARNING: if your computer has fewer than ",nruns[2]," processor cores, this will do more harm than good\n",sep="")
	} else {
		parallel<-FALSE
	}

	# pre-compute edit cost matrix / all dev seq mat 
	# and pass to parallel runs otherwise, it will be RE-computed for each run
	if(pgi.tree$inf.params$heuristic == "exhaustive") {
		if(verbosity > 1) cat("Precomputing all developmental sequence matrix...\n")
		pgi.tree$inf.params$all.dev.seq.mat<-unique(rbind(create.all.dev.seq.mat(pgi.tree$seq.len,pgi.tree$inf.params$simultaneity,verbosity=verbosity) ,extract.seqs.from.tree(pgi.tree$tree,pgi.tree$seq.len,c())))
	}
	if(pgi.tree$inf.params$heuristic == "sw04") {
		if(verbosity > 1) cat("Precomputing edit cost matrix...\n")
		t.edmat<-create.edit.cost.matrix(pgi.tree,verbosity=verbosity)
		pgi.tree$inf.params$edit.cost.matrix<-t.edmat[[2]]
		pgi.tree$inf.params$all.dev.seq.mat<-t.edmat[[1]]
	}
	if(verbosity == 0) pb<-txtProgressBar(min=0,max=total.runs,style=3)
	if(all(nruns == c(1,1))) {
		# a single run, the default
		if(!(pgi.tree$inf.params$resume && (pgi.tree$type == "filled"))) {
			pgi.tree<-pgi.inference(pgi.tree, verbosity=verbosity)
		} 
		if(pgi.tree$inf.params$resume) {
			partial.pgi.tree<-pgi.tree
			save(partial.pgi.tree,file=pgi.tree$inf.params$resume.temp.file)
			rm("partial.pgi.tree")
		}
		pgi.con.trees<-list(pgi.pseudo.consensus(pgi.tree,verbosity=verbosity))
		if(pgi.tree$inf.params$resume) unlink(pgi.tree$inf.params$resume.temp.file)
	} else {
		pgi.con.trees<-list()
		nseries<-nruns[1]
		num.par<-nruns[2]
		t.filled.trees<-list()
		for(i in 1:nseries) {
			if(parallel) {
				require(foreach)
				require(doMC)
				registerDoMC(cores=num.par)	
				init<-foreach(j=1:num.par,.export=c("verbosity","num.par")) %dopar% (paste("Initializing parallel processing, core",j,"of",num.par,"\n"))
				if(verbosity > 0) print(init)
				if(verbosity < 2)	verb<-c(rep(0,times=num.par-1),verbosity) else verb<-rep(verbosity,times=num.par)
				cur.trees<-seq((length(t.filled.trees)+1),(length(t.filled.trees)+num.par))
				t.filled.trees<-c(t.filled.trees,foreach(j=1:num.par,.export=c("pgi.tree","verb","num.par"),.packages=c("e1071","ape","pgi2")) %dopar% (pgi.inference(pgi.tree=pgi.tree, verbosity=verb[j])))
				pgi.con.trees<-c(pgi.con.trees,foreach(j=1:num.par,.export=c("pgi.tree","verb","num.par"),.packages=c("e1071","ape","pgi2")) %dopar% (pgi.pseudo.consensus(pgi.tree=t.filled.trees[[cur.trees[j]]],verbosity=verb[j])))
				if(verbosity == 0) setTxtProgressBar(pb,getTxtProgressBar(pb)+num.par)
			} else {
				pgi.tree<-pgi.inference(pgi.tree, verbosity=verbosity)
				pgi.con.trees<-c(pgi.con.trees,list(pgi.pseudo.consensus(pgi.tree,verbosity=verbosity)))
				if(verbosity == 0) setTxtProgressBar(pb,getTxtProgressBar(pb)+1)
			}
			if(verbosity > 0) cat("Completed cycle",i,"out of",nruns[1]," cycles in series\n")
		}
	}
	class(pgi.con.trees)<-"multi.pgi.tree"
	pgi.tree.length.dist<-sapply(pgi.con.trees,function(x){return(x$tree.length)})	
	if(verbosity > 0) {
		cat("\nResults\n")
		cat("Minimum tree length = ",min(pgi.tree.length.dist),sep="")
		if(sum(nruns) > 1) {
			cat(", mean tree length across replicates = ",mean(pgi.tree.length.dist),"\n",sep="")
		} else {
			cat(" based on a single run\n")
		}
	}
	if(verbosity == 0) close(pb)
	if(interactive) pgi.interactive.consensus(pgi.con.trees)
	cat("Executed cleanly\n")
	if(supercon && (length(pgi.tree.length.dist) > 1)) {
		if(verbosity > 1) cat("Calculting super-pseudo-consensus tree using a tolerance of ",(if(is.na(supercon.tol)) "NA - use all trees" else supercon.tol),"\n")
		return(pgi.supercon(pgi.con.trees,tol=supercon.tol,verbosity=verbosity))
	} else {
		if(length(pgi.tree.length.dist) > 1) {
			return(pgi.con.trees)
		} else return(pgi.con.trees[[1]])
	}
}

