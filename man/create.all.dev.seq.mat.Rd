\name{create.all.dev.seq.mat}
\alias{create.all.dev.seq.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO CREATE A MATRIX OF ALL POSSIBLE DEVELOPMENTAL SEQUENCES FOR A GIVEN NUMBER OF ELEMENTS
}
\description{
	This function returns a matrix of all possible ranked developmental sequenes (rows) for a given sequence length. 
}
\usage{
create.all.dev.seq.mat(seq.len, simul, verbosity = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq.len}{
		Number of events (sequence length)
}
  \item{simul}{
		TRUE/FALSE: should sequences including simultaneity be included?
}
  \item{verbosity}{
		0 - no output, 1 - progress bar, 2+ diagnostic output.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A matrix of ranked developmental sequences (rows)
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
	\code{\link{create.edit.cost.matrix}}
}
\examples{
	all.dev.seq.mat<-create.all.dev.seq.mat(seq.len=4,simul=TRUE,verbosity=1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
