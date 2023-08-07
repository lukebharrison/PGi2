\name{pgi.con.tree.length}
\alias{pgi.con.tree.length}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO DETERMINE THE TREE LENGTH OF A PGI CONSENSUS TREE
}
\description{
 WARNING: only results from the "simple" pseudoconsensus will make sense with this function. Returns the tree length in terms of sequence heterochronies of a given PGi tree datastructure.
}
\usage{
pgi.con.tree.length(con.tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con.tree}{
		A PGi pseudoconsesnus datastructure (e.g. pgi.tree$consensus)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 Tree length, in number of sequence heterochronies
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
	
	data(velhagen.simple.con)
	pgi.con.tree.length(velhagen.simple.con$consensus)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
