\name{pgi.assemble.pgi.tree}
\alias{pgi.assemble.pgi.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	UTILITY FUNCTION TO CONVERT RAW PGI TREE DATA STRUCTURES (OLD) TO THE PGI DATA OBJECT
}
\description{
	Taking an old PGi tree data structure, creates the PGi data object around it.
}
\usage{
pgi.assemble.pgi.tree(old.pgi.tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{old.pgi.tree}{
		An old PGi tree data structure (or pgi.tree$tree, for fun)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A PGi tree data object of class "pgi.tree"
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

	data(velhagen)
	int.tree<-velhagen$tree
	velhagen.new<-pgi.assemble.tree(int.tree)
	# velhagen should now be identical to velhagen.new	


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
