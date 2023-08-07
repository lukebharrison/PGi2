\name{seq.to.rds}
\alias{seq.to.rds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FUNCTION TO CONVERT A DEVELOPMENTAL SEQUENCE INTO A RANKED DEVELOPMENTAL SEQUENCE
}
\description{
Developmental sequences are stored as lists of elements, with each element a vector of one or more events. If the number of events per list element is greater than one, this indicates simulateity. This function converts the list-based developmental sequence structure into a ranked developmental sequence, which consists of a vector of length = number elements and for each element, its rank order.
}
\usage{
seq.to.rds(dev_seq, seq_len)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dev_seq}{
		A developmental sequence in PGi format - a list of integer vectors (length >= 1).
}
  \item{seq_len}{
		The total number of elements present (ranked developmental sequence length)
}
}
\details{
	Note: PGi stores developmental sequence as lists and ranked developmental seqeunces as vectors.
}
\value{
	The developmental sequence in ranked form.
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
	\code{\link{rds.to.seq}}
}
\examples{

	dev.seq<-list(1,c(2,3),4)
	rds.seq<-seq.to.rds(dev.seq,4)
	


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
