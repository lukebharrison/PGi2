\name{compress.rds}
\alias{compress.rds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO "COMPRESS" A RANKED DEVELOPMENTAL SEQUENCE
}
\description{
	This function "compresses" a ranked developmental sequence to the fewest possible ranks while preserving rank-order.
}
\usage{
compress.rds(rds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds}{
	A ranked developmental sequence of length = number of elements with corresponding ranks
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A "compressed" ranked developmental sequence using the smallest number of ranks possible while retaining all rank-order information.
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

	rds<-matrix(c(1,7,2,3,4),1,5)
	rds<-compress.rds(rds)
	

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
