\name{create.edit.cost.matrix}
\alias{create.edit.cost.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO CREATE AN EDIT COST MATRIX FOR SEARCH-BASED OPTIMIZATION
}
\description{
	This function generates an edit-cost matrix to perform search based optimization: the pair-wise edit cost between all possible sequences is calculated based on the chosen edit cost function in the inf.params data structure (or simple list(edit.cost.func="...",simultaneity=TRUE)).
}
\usage{
create.edit.cost.matrix(pgi.tree, inf.params = NULL, verbosity = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pgi.tree}{
		A pgi.tree data object (to ensure the tip seqeunce are included, even if simulataneous sequences are not)
}
  \item{inf.params}{
		The inference parameters data structure (a named list). Minimally, provide list(edit.cost.func=X,simultaneity=Y).
}
  \item{verbosity}{
		0 - no output, 1 - progress bar, 2+ diagnostic output.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
	\code{\link{create.all.dev.seq.mat}}
}
\examples{



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
