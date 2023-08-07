pgi.interactive.params <-
function(pgi.tree) {
	cat("\nSet up the parameters for the analysis\n\n")
	cat("Number of independent PGi computations to be run. It is recomended to never run less than ~4 computations because PGi is a heuristic method\n")
	cat("Note: If you will be running computations in parallel on a multi-processor/core computer, please choose a multiple of the number of cores to be used\n")
	nruns<-as.integer(readline(prompt="Number of runs? > "))
	cat("Please choose the number of runs to execute simultaneously in parallel (Warning: do not exceed the number of cores in this computer, it is counter-productive)\n")
	n.par<-as.integer(readline(prompt="Number of run[s] to execute simultaneously? > "))
	nruns<-c(round(nruns/n.par),n.par)
	cat("PGi will execute ",nruns[1]*nruns[2]," run[s], in ",nruns[1]," iterations of ",nruns[2]," parallel run[s].\n\n",sep="")
	pgi.tree$nruns<-nruns
	# heuristic parameters
	cat("\nWhich heuristic to use?.\n- 'pgi' for standard PGi algorithm. (use this)\n- 'oldgenetic' for a very similar heuristic that does not stop inference at a node if it detects no new sequences are being discovered.\n- 'exhaustive' for a non-heuristic method that samples *ALL* possible ancestral sequences at each node (very, very slow - only use in very small data sets).\n- 'sw04' for Search-based optimization, as used in Schulmeister and Wheeler (2004), also select 'event-pair' for the edit-cost function below.\n\n")
	heuristic<-readline(prompt="Heuristic? > ")
	if(any(heuristic == c("pgi","oldgenetic"))) {
		cat("Parameters for the PGi Simplified genetic algorithm\n Note: For each of the following variables, an increase in the parameter will increase the accuracy of the analysis but also increase running time")
		cat("Number of cycles of 'selection' per node: best to choose an integer of at least 2*(sequence length)), empirical cycles = 100 for a sequence length of 24 has worked\n") 
		cycles<-as.integer(readline(prompt="Number of cycles of 'selection' per node? > "))
		cat("Number of replicate perturbed sequences per cycle of 'selection' per node (e.g. variation): best to choose an integer of at least 2*(sequence length))\n") 
		replicates<-as.integer(readline(prompt="Number of replicates per cycle of 'selection' per node? > "))
		if(heuristic == "pgi") {
			cat("Maximum number of ancestral developmental sequence to be retained at each node: also recomended at least 2*(sequence length). A value of 0 retains all sequences\n")
		ret.anc.seq<-as.integer(readline(prompt="Retained number of ancestral sequences per node? > "))
		}
	}
	#edit cost function
	cat("\nWhich edit-cost function to use?\n- 'parsimov' for the Parsimov algorithm (Jeffery et al 2005) with a greedy heuristic (use this)\n- 'exhaustive-parismov' for the Parsimov algorithm but without the greedy heuristic\n- 'event-pair' for event-pair distance (e.g. Jeffery et al. 2002; Schulmeister and Wheeler, 2004)\n")
	edit.cost.func<-readline(prompt="Edit cost function? > ")
	pgi.tree$inf.params<-list(replicates=replicates,cycles=cycles,ret.anc.seq=ret.anc.seq,simultaneity=TRUE,heuristic=heuristic,edit.cost.func=edit.cost.func)
	# consensus methodologies
	cat("\nWhich Pseudocosensus methodology to apply to the results?\n- 'exhaustive' make a consensus of all possible solutions of equal score, at each node in the phylogeny (use this if possible).\n- 'semi-exhaustive' make a consensus of all possible solutions of equal score, at each node, but with a maximum number of evaluated solutions (use this if exhaustive stalls).\n- 'simple' consensus of all solutions of equal score at the root note only (very fast, but must be used with *many* independent runs (use only if the first two are too slow).\n")
	con.type<-readline(prompt="Consensus method? > ")
	if(con.type == "semi-exhaustive") {
		cat("\nPlease select a maximum number (per node) of solutions to evaluate before stopping, this will be split among the equally good solutions at the root. Choose a large value, say 2000-5000?\n")
		semi.ex.con.max.n<-as.integer(readline(prompt="Max number of solutions? > "))
	} else { 
		semi.ex.con.max.n<-0
	}
	pgi.tree$con.params<-list(con.type=con.type,semi.ex.con.max.n=semi.ex.con.max.n,edit.cost.func=edit.cost.func)

	cat("\n\nVerbosity?\n-'0': a sigle progress bar for all runs.\n-'1': seperate progress bars for each run + consensus.\n-'2' or more: increasing amount of diagnostic output, recomended only for troubleshooting.\n")
	verbosity<-as.integer(readline(prompt="Verbosity? > "))
	pgi.tree$verbosity <- verbosity
	if(match.arg(readline(prompt="\nDisplay the PGi data object to verify topology + parameters? > "),c("yes","no")) == "yes") {
		cat("\n")		
		print(summary(pgi.tree,print.all.seqs=TRUE))
		cat("\n")
		if(match.arg(readline(prompt="\nIs this correct? > "),c("yes","no")) == "no") stop("Incorrect parameters, please run pgi.interactive() again\n")
	}
	return(pgi.tree)
}

