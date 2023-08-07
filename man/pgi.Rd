\name{pgi}
\alias{pgi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	MAIN PGI CONTROL FUNCTION
}
\description{
The pgi() function handles execution of the PGi analysis and pseudoconsensus. This function can run interactively (interactive=TRUE) or non-interactively (interactive=FALSE). PGi can be executed in series or in parallel (WARNING: Not working on all computers) using the foreach() function from the foreach and doMC libraries. 
}
\usage{
pgi(pgi.tree, interactive = FALSE, nruns = NULL, verbosity = 0, replicates = 100, cycles = 100, ret.anc.seq = 100, simul = TRUE, heuristic = "pgi", edit.cost.func = "parsimov", con.type = "exhaustive", semi.ex.con.max.n = 0, supercon = FALSE, supercon.tol = NA, overwrite = FALSE, inf.params = NULL, con.params = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pgi.tree}{
		PGi data object of class "pgi.tree" which contains the topology and developmental sequences - can be left blank with interactive=TRUE: you will prompted to load from disk.
}
  \item{interactive}{
		TRUE/FALSE: executive interactively? If TRUE, you will be prompted 
}
  \item{nruns}{
		Vector of 2 integers. The first integer specifies the number of cycles of execution in series. The second interger specifies the number of executions in parallel. The total number of runs is nruns[1]*nruns[2].
}
  \item{verbosity}{
		Vebosity: 0: overall progress bar for all runs, 1: progress bar for each independent run, 2: progress bar for each node of each independent run, 3+ diagnostic output
}
  \item{replicates}{
		See inf.params and details, below
}
  \item{cycles}{
		See inf.params and details, below
}
  \item{ret.anc.seq}{
		See inf.params and details, below
}
  \item{simul}{
		See inf.params and details, below
}
  \item{heuristic}{
		See inf.params and details, below
}
  \item{edit.cost.func}{
		See inf.params and details, below
}
  \item{con.type}{
		See con.params and details, below
}
  \item{semi.ex.con.max.n}{
		See con.params and details, below
}
  \item{supercon}{
		TRUE/FALSE: automatically generate a superconsensus tree of the pseudoconsensus trees from each run? see ?pgi.supercon for more information.
}
  \item{supercon.tol}{
		Float [0-1]: When generating the superconensus tree, normally all trees, regardless of tree length are included. If specified, this is the fraction of deviation from the minimum tree length that will be included: 0 - all pseudoconsensus trees includes, 1 - only independent runs with the minimum observed tree length are included.
}
  \item{overwrite}{
		TRUE/FALSE: if TRUE, if pgi.tree object that already contains a consensus is provided, the $filled and $consensus slots will be overwritten, otherwise an error will be thrown.
}
  \item{inf.params}{
		If inf.params is not provided in the pgi.tree object, provide here, or override the the $inf.params provided in the pgi.tree with an explicit specification.
}
  \item{con.params}{
		If con.params is not provided in the pgi.tree object, provide here, or override the the $inf.params provided in the pgi.tree with an explicit specification.
}
}
\details{
The inference proceedure is controlled by a set of argument normally contained in the inf.params data structure. This is a named list of parameters:
	inf.params$heuristic 				- string:		the (heuristic) method for inference. choice of "pgi", "oldgenetic", "exhaustive" or "sw04" (details below).
	inf.params$cycles						- integer:	number of cycles of selection for the genetic heuristic at each node (for heuristic = "pgi" || "oldgenetic").
	inf.params$replicates				- integer: 	number of replicate hypothetical ancestral sequences to generate per cycle of selection per node (for heuristic = "pgi" || "oldgenetic").
	inf.params$ret.anc.seq			- integer:	number of unique hypothetical ancestral sequences to retain at each node.
	inf.params$edit.cost.func		- string:		edit cost function to apply to evaluate hypothetical anctreal sequences. choice of "parsimov", "exhaustive-parsimov" or "event-pair" see ?edit.cost for details.
	inf.params$simultaneity			- TRUE/FALSE:	If FALSE, consider only fully resolved ancestral sequences. If TRUE, consider both fully resolved and sequences with simultaneity.
	inf.params$all.dev.seq.mat	-	matrix:		for heuristic == "exhaustive" || "sw04", the matrix of all possible developmental sequences for a given length (will be created if missing).
	inf.params$edit.cost.matrix	- matrix: 	for heuristic == "sw04", the edit-cost matrix of edit costs between all possible sequences (will be created if missing).

	(Heuristic) Method Types:
	"pgi" the algorithm described in Harrison and Larsson(2008). Uses a simplified genetic algorithm to refine hypothetical ancestral sequences at each node. Applies an additional check for convergence to save computation: if addition compution no longer improves the hypothetical ancestral sequence, execution at this node is truncated.
	"oldgenetic": same as "pgi" but without the addition check: execution continues for the full amounts specified by cycles and replicates.
	"exhaustive": try each and every possible ancestral sequence at each node. Slow but non-heuristic.
	"sw04":	search-based optimization as described in Schulmeister and Wheeler (2004). Precomputes pair-wise edit costs between all possible sequences and then transverses the tree. Faster and smarter non-heuristic compared with the "exhaustive" approach.

	Pseudoconsensus parameters are specific in the con.params datastrucutre, a named list of parameters. This is either provided in the pgi.tree data object or explicitely as an argument. Data structure list members:
	con.params$con.type: 					String: the pseudoconsensus method to use. Choice of "exhaustive", "semi-exhaustive" or "simple" (see descriptions in notes, below)
	con.params$semi.ex.con.max.n 	Integer: In the semi-exhaustive pseudoconsensus, the maximum (per branch) transversals. A number ~ 2000 is reasonable and prevents excessive computation.
	con.params$edit.cost.func			String: The edit cost function to compute sequence heterochronies. Must be one of "parsimov" or "exhaustive-parsimov" even if the ancestral sequences were inferred under event-pair distance.

	PGi Pseudoconsensus Methods:
	"exhaustive": exhaustive pseudoconsensus - Following traceback pointers from each "equally parsimonious" hypothetical ancestral developmental sequence at the root of the tree, the pseudoconsensus method, at every node, examines every possible combination of ancestral and descendant sequences (from those recorded during inference) and records any other possible solutions for that subtree that have an equal score. In this fashion, the pseudoconesus method explores a degree of solution space to determine how often specific sequence heterochronies are recovered.
	"semi-exhaustive": semi-exhaustive pseudoconsensus - nearly identical to the "exhaustive" method, but includes a cut-off for the number of transversal per branch. This is necessary, as the exhaustive method will revisit the same branches many times, espeically near the tips and become for all intents and purposes stuck.
	"simple": in contrast to the exhaustive and semi-exhaustive methods, the simple pseudoconsensus method simply follows the traceback pointers from each "equally parsimonious" solution at the root to the tips but does not explore alternative and equally scoring solutions for each subtree. This should only be used in conjunction with multiple indenpendent runs and integrated with a superconsensus method.


}
\value{
	A data object of type "multi.pgi.tree" consisting of list of pseudoconsensus trees generated by the pgi() function. If supercon==TRUE, then a single object of class "pgi.tree" is returned in the superconsensus tree.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pgi.inference}},\code{\link{pgi.pseudoconsensus}},\code{\link{pgi.supercon}}
}
\examples{

	# fully interactive
	pgi(interactive=TRUE)

	# interactive
	data(velhagen)
	con.trees<-pgi(velhagen,interactive=TRUE,supercon=FALSE)

	# non-interactive
	data(sanchez.villagra.2002)
	con.trees<-pgi(sanchez.villagra.2002,nruns=c(2,1),inf.params=list(heuristic="pgi",cycles=100,replicates=100,ret.anc.seq=100,edit.cost.func="parsimov",simultaneity=TRUE),con.params=list(con.type="semi.exhaustive",semi.ex.con.max.n=5000,edit.cost.func="parsimov"),verbosity=1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
