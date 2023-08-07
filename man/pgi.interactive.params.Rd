\name{pgi.interactive.params}
\alias{pgi.interactive.params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO INTERACTIVELY CHOOSE PARAMETERS FOR PGI EXECUTION
}
\description{
	This function is called on a pgi.tree object in order to store execution parameters for the pgi analysis. Using the console the user is invited to paramterize the inference and pseudoconsensus stages of the pgi analysis. These are store in special fields in the pgi.tree object which make it ready to be executed using pgi().
}
\usage{
pgi.interactive.params(pgi.tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pgi.tree}{
		A pgi.tree data object containing topology and developmental sequences. Usually this is created using pgi.read.nexus() or pgi.assemble.tree()
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	An identical pgi.tree data object, but with $con.params, $inf.params, $nruns and $verbosity fields updated or created.
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
\code{\link{pgi}},\code{\link{pgi.inference}},\code{\link{pgi.pseudoconsensus}}
}
\examples{
	

	data(velhagen)
	velhagen<-pgi.interactive.params(velhagen)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
