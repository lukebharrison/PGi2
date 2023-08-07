\name{rds.to.seq}
\alias{rds.to.seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FUNCTION TO CONVERT A RANKED DEVELOPMENTAL SEQUENCE INTO A DEVELOPMENTAL SEQUENCE
}
\description{
Ranked developmental sequences are represented as integer vectors with a lenght of the number of elements. For each position i, element i's ranked is stored. In the case of missing/unscorable elements, -2 is recorded. This function coverts a ranked developmental sequence into a development sequence. Developmental sequences in PGi are stored as lists of elements, with each element a vector of one or more events. If the number of events per list element is greater than one, this indicates simulateity.
}
\usage{
rds.to.seq(rds_dev_seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds_dev_seq}{
		A ranked developmental sequence of length = number of elements with corresponding ranks
}
}
\details{
	Note: Missing values are recorded with -2
}
\value{
	The developmental sequence in PGi's developmental sequence format (list of elements)
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
	\code{\link{seq.to.rds}},\code{\link{print.rds}},\code{\link{compress.rds}}
}
\examples{

	rds.seq<-matrix(c(4,3,3,2,1),1,5)
	dev.seq<-rds.to.seq(rds.seq)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
