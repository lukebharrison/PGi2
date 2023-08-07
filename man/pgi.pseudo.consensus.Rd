\name{pgi.pseudo.consensus}
\alias{pgi.pseudo.consensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO CALCULATE THE PGI PSEUDOCONSENSUS OF A FILLED PGI DATA OBJECT
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pgi.pseudo.consensus(pgi.tree, verbosity = 0, con.params = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pgi.tree}{
		A pgi.tree data object containing minimally a $filled.tree slot for pseudoconsensus and usually, a $con.params slot with pseudoconsensus parameters.
}
  \item{verbosity}{
		Verbosity 0: no outputs, 1: overall progress bar, 2: per node progress bar, 3+ diagnostic output
}
  \item{con.params}{
		The pseudoconsensus paramters data structure. Normally, $con.params is provided in the pgi.tree object. However, if this is not the case or you wish to override the values stored in the pgi.tree data object, you can specify alternative paramters. See details, below for a description of the con.params data structure.
}
}
\details{
	Pseudoconsensus parameters are specific in the con.params datastrucutre, a named list of parameters. This is eithe provided in the pgi.tree data object or explicitely as an argument. Data structure list members:
	con.params$con.type: 					String: the pseudoconsensus method to use. Choice of "exhaustive", "semi-exhaustive" or "simple" (see descriptions in notes, below)
	con.params$semi.ex.con.max.n 	Integer: In the semi-exhaustive pseudoconsensus, the maximum (per branch) transversals. A number ~ 2000 is reasonable and prevents excessive computation.
	con.params$edit.cost.func			String: The edit cost function to compute sequence heterochronies. Must be one of "parsimov" or "exhaustive-parsimov" even if the ancestral sequences were inferred under event-pair distance.
}	
\value{
	An object of class "pgi.tree" identical to the input object, but updated with: $con.params, if explicitely provided and $consensus which holds the pseudoconsensus data structure. Additionally, the $type field is updated to "consensus"
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
	PGi Pseudoconsensus Methods:
	"exhaustive": exhaustive pseudoconsensus - Following traceback pointers from each "equally parsimonious" hypothetical ancestral developmental sequence at the root of the tree, the pseudoconsensus method, at every node, examines every possible combination of ancestral and descendant sequences (from those recorded during inference) and records any other possible solutions for that subtree that have an equal score. In this fashion, the pseudoconesus method explores a degree of solution space to determine how often specific sequence heterochronies are recovered.
	"semi-exhaustive": semi-exhaustive pseudoconsensus - nearly identical to the "exhaustive" method, but includes a cut-off for the number of transversal per branch. This is necessary, as the exhaustive method will revisit the same branches many times, espeically near the tips and become for all intents and purposes stuck.
	"simple": in contrast to the exhaustive and semi-exhaustive methods, the simple pseudoconsensus method simply follows the traceback pointers from each "equally parsimonious" solution at the root to the tips but does not explore alternative and equally scoring solutions for each subtree. This should only be used in conjunction with multiple indenpendent runs and integrated with a superconsensus method.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pgi}},\code{\link{pgi.supercon}}
}
\examples{

	data(velhagen.filledtree)
	con.trees<-pgi.pseudoconsensus(velhagen.filledtree,con.params=list(con.type="exhaustive",edit.cost.func="parsimov"),vebosity=1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
