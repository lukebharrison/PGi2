\name{pgi.inference}
\alias{pgi.inference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	PERFORM THE PGI INFERENCE STEP TO CHOOSE RECONSTRUCT ANCESTRAL DEVELOPMENTAL SEQUENCES
}
\description{
	Note: Please use the pgi() main function to call this function unless you know what you are doing.
	pgi.inferece() is the main control function for the first step of the PGi analysis. This function controls the parameters for the analysis to fill the empty PGi data object with hypothetical ancestral sequences. Inference execution parameters are normally provided by the inclusion of a inf.params data structure (list of named parameters) within the pgi.tree object. Alternative, these can be overridden by explicitely provided an inf.param data structure as an argument (details below). 
}
\usage{
pgi.inference(pgi.tree, verbosity = 1, inf.params = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pgi.tree}{
		PGi data object. If the data object already contains a filled or pseudoconsensus tree, it will be overriden. 
}
  \item{verbosity}{
		Verbosity: 0 : no messages, 1 : overall progress bar, 2: per node progress bars, 3+ diagnostic output.
}
  \item{inf.params}{
		If the PGi data object does not include an inference parameters (inf.params) data structure, one must be specific. Providing this argument will override any inf.params structure in the PGi data object. see details below.
}
}
\details{
	The inference proceedure is controlled by a set of argument normally contained in the inf.params data structure. This is a named list of parameters:
	inf.params$heuristic 				- string:		the (heuristic) method for inference. choice of "pgi", "oldgenetic", "exhaustive" or "sw04" (details below).
	inf.params$cycles						- integer:	number of cycles of selection for the genetic heuristic at each node (for heuristic = "pgi" || "oldgenetic").
	inf.params$replicates				- integer: 	number of replicate hypothetical ancestral sequences to generate per cycle of selection per node (for heuristic = "pgi" || "oldgenetic").
	inf.params$ret.anc.seq			- integer:	number of unique hypothetical ancestral sequences to retain at each node.
	inf.params$edit.cost.func		- string:		edit cost function to apply to evaluate hypothetical anctreal sequences. choice of "parsimov", "exhaustive-parsimov" or "event-paur" see ?edit.cost for details.
	inf.params$simultaneity			- TRUE/FALSE:	If FALSE, consider only fully resolved ancestral sequences. If TRUE, consider both fully resolved and sequences with simultaneity.
	inf.params$all.dev.seq.mat	-	matrix:		for heuristic == "exhaustive" || "sw04", the matrix of all possible developmental sequences for a given length (will be created if missing).
	inf.params$edit.cost.matrix	- matrix: 	for heuristic == "sw04", the edit-cost matrix of edit costs between all possible sequences (will be created if missing).
}
\value{
	A "pgi.tree" class object identical to the input pgi.tree but including: 
	- updated inf.params data structure if explicitely specified
	- type now == "filled"
	- a pgi.tree$filled.tree data structure that contains the hypothetical ancestral sequences and their scores, ready for pseudoconsensus analysis.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
	(Heuristic) Method Types:
	"pgi" the algorithm described in Harrison and Larsson(2008). Uses a simplified genetic algorithm to refine hypothetical ancestral sequences at each node. Applies an additional check for convergence to save computation: if addition compution no longer improves the hypothetical ancestral sequence, execution at this node is truncated.
	"oldgenetic": same as "pgi" but without the addition check: execution continues for the full amounts specified by cycles and replicates.
	"exhaustive": try each and every possible ancestral sequence at each node. Slow but non-heuristic.
	"sw04":	search-based optimization as described in Schulmeister and Wheeler (2004). Precomputes pair-wise edit costs between all possible sequences and then transverses the tree. Faster and smarter non-heuristic comapred with the "exhaustive" approach.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pgi}}
}
\examples{

	data(velhagen)
	inf.params<-list(heuristic="pgi",cycles=50,replicates=50,ret.anc.seq=50,edit.cost.func="parsimov",simultaneity=TRUE)
	pgi.inference(velhagen,inf.params=inf.params,verbosity=1)
	


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
