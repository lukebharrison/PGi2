pgi.pseudo.consensus <-
function(pgi.tree,verbosity=0,con.params=NULL) {
	if(pgi.tree$type != "filled") {
		if(pgi.tree$type == "empty") stop("Can't make a pseudo consensus of an empty PGi tree, run infer() first.")
		if(pgi.tree$type == "consensus") {
			cat("Warning: overwriting existing consensus tree, interrupt now to abort (2 seconds)\n")
			Sys.sleep(2)
		}
	}		
	if(!is.null(con.params)) {
		pgi.tree$con.params<-con.params
	} else {
		if(is.null(pgi.tree$con.params)) stop("No consensus parameters stored in the pgi.tree object nor provided, stopping..")
	}
	if(pgi.tree$con.params$con.type == "exhaustive") {
		pgi.tree$con.params$semi.ex.con.max.n<-0
		if(verbosity > 0) cat("Selected exhaustive pseudoconsensus. WARNING: on some data sets, this consensus method can stall.\n")
		return(.pgi.exhaustive.consensus(pgi.tree,verbosity=verbosity))
	} else if(pgi.tree$con.params$con.type == "semi-exhaustive") {
		if((pgi.tree$con.params$semi.ex.con.max.n)==0) stop("Please provide a per branch maximum transversal for the semi-exhaustive pseudoconsensus (e.g. semi.ex.con.max.n = 5000)")
		if(verbosity > 0) cat("Selected semi-exhaustive pseudoconsensus.\n")
		return(.pgi.exhaustive.consensus(pgi.tree,verbosity=verbosity))
	} else 	if(pgi.tree$con.params$con.type == "simple") {
		if(verbosity > 0) cat("Selected simple pseudoconsensus. WARNING: exhaustive or semi-exhaustive is prefered, if computable\n")
		return(.pgi.simple.con(pgi.tree,verbosity=verbosity))
	} else stop("Invalid consensus type: choose among 'exhaustive', 'semi-exhaustive' or 'simple'")
}

