\name{pgi.interactive.consensus}
\alias{pgi.interactive.consensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO INTERACTIVELY DESCRIBE THE RESULTS OF THE PGI ANALYSIS AND PRODUCE PLOTS
}
\description{
	Note: best used via the pgi(interactive=T) function.
	This function opens an interactive session to analyze and plot the results of a PGi analysis. Based on an object of class "multi.pgi.phylo", this function interactively asks whether the user would like to generate and plot a superconsensus of pseudoconsensus trees.
}
\usage{
pgi.interactive.consensus(con.trees, verbosity = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con.trees}{
		An object of class "multi.pgi.tree" which is a list of pgi.tree objects, in this case, they must contain pseudoconsensus trees in the $consensus slot (derived from the pgi.pseudoconsensus.function)
}
  \item{verbosity}{
		Leave at default, for internal use.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	No explict return() but this function will assign to variables of the user's choosing during execution.
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
\code{\link{pgi.pseudoconsensus}},\code{\link{pgi}}
}
\examples{

	data(velhagen.pseudocontrees)
	pgi.interactive.consensus(velhagen.pseudocontrees)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
