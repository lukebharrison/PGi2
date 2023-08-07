pgi.interactive.consensus <-
function(con.trees,verbosity=1) {
	cat("Provide an object name (file will also be saved to HD) for the PGi data object containing the pseudoconsensus trees")
	fn<-readline(prompt="Object name? > ")
	assign(fn,con.trees,envir = .GlobalEnv)
	save(con.trees,file=paste(fn,".rda",sep=""))

	plt<-match.arg(readline(prompt="\nPlot the [first/only] run's pseudoconsensus tree? Type 'yes' to plot or no to continue > "),c("yes","no"))
	if(plt == "yes") {
		if(class(con.trees) == "multi.pgi.tree") {
			plot(con.trees[[1]]) 
		} else {
			plot(con.trees)
		}
	}
	cat("\n")		
	if(class(con.trees) == "multi.pgi.tree") {
		scon<-match.arg(readline(prompt="Generate a super consensus tree? >"),c("yes","no"))
		if(scon == "yes") {
			supercon.tree<-pgi.supercon(con.trees,verbosity=verbosity)
			cat("Provide an object name (file will also be saved to HD) for the PGi data object containing the SUPER pseudoconsensus trees\n")
			fn<-readline(prompt="Object name? > ")
			assign(fn,supercon.tree,envir = .GlobalEnv)
			save(supercon.tree,file=paste(fn,".rda",sep=""))
			plt<-match.arg(readline(prompt="\nPlot the super pseudoconsensus trees? [no]> "),c("yes","no",""))
			if(plt == "yes") {
				plot(supercon.tree)
			}
		} else {
			cat("Superconsensus not generated (use pgi.supercon() on the ",fn," object later if you want.\nQuiting.\n",sep="")
		}
	} else {
			cat("Only a single consensus tree here.\n")
	}
}

