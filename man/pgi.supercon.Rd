\name{pgi.supercon}
\alias{pgi.supercon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO CALCULATE A SUPERCONSENSUS OF PGI PSEUDOCONSENSUS TREES
}
\description{
	If mulitple independent runs of PGi have been calcualted, this function allows the user to combine the pseudoconsensus trees into a single super consensus tree. Basically, this is an average of the sequence heterochronies and ancestral sequences recovered in each pseudoconsensus analysis. The results is returned as a "pgi.tree" object of type "superconesnsus"
}
\usage{
pgi.supercon(con_trees, tol = 0, verbosity = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con_trees}{
	A list of pgi.tree data objects. This should be of class "multi.pgi.tree".
}
  \item{tol}{
	Tolerance: the maximum deviation in terms of sequence heterochronies from the shortest input pseudoconsensus to be used for the superconsensus. -1: use all.
}
  \item{verbosity}{
	Verbosity: 0 - no output to console, 1 - simple messages, 2+ diagnostic output.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A "pgi.tree" data object with the super conesnsus stored in the consensus slot and type = "superconsensus".
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
