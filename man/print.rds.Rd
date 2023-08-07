\name{print.rds}
\alias{print.rds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO FORMAT A RANKED DEVELOPMENTAL SEQUENCE FOR PRINTING
}
\description{
Function to format a ranked development sequence (a vector of integers) for printing. Retunrns the sequence using square brackets to indicate simultaneity and commas to seperate ranks.
}
\usage{
print.rds(rds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds}{
		A ranked developmental sequence in PGi format (vector of integers of ranks for elements)
}
}
\details{
	Note: PGi represents ranked developmental sequences as vectors of length = number of elements, where the i-th integer represents the rank of the i-th element.
}
\value{
	A string that consists of the nicely formatted developmental sequence.
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
	\code{\link{seq.to.rds}},\code{\link{rds.to.seq}},\code{\link{compress.rds}}
}
\examples{

	rds.seq<-matrix(c(1,3,3,2,1),1,5)
	print(print.rds(rds.seq))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
