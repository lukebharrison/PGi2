\name{edit.cost}
\alias{edit.cost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO COMPUTE EDIT COST BETWEEN TWO RANKED DEVELOPMENTAL SEQUENCES
}
\description{
	This function computes the edit cost between two developmental sequences. In this implementation, between an ancestor and a descendant sequence. 
}
\usage{
edit.cost(anc_dev_seq, der_dev_seq, output, def_penalty, method, seq_len)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{anc_dev_seq}{
		The Ancestral ranked developmental sequence (integer vector)
}
  \item{der_dev_seq}{
		The descedant ranked developmental sequence (integer vector)
}
  \item{output}{
		This flag determines how the output is returned (only for parsimov or exhaustive-parsimov methods):
		- 0:	Number of sequence heterochronies + missing elements "penalty" (see notes below)
		- 1:	Number of sequence heterochronies
		- 2:	Vector of sequence heterochronies, returned with element number, negative values are deccelerations, positive values accelerations
		- 3:  A list of two elements: [[1]]: as #2, the vector of infered sequence heterochronies and [[2]] a vector identifying elements "lost" or "created" along the branch
}			
  \item{def_penalty}{
		A tuning parameter for missing data (leave unspecified unless you have a compelling reason to change)
}
  \item{method}{
		The edit cost function to use: 
		-	"parsimov":							the Parsimov function as described in Jeffery et al (2005), using a greedy heuristic
		- "exhaustive-parsimov": 	the Parsimov function as described in Jeffery et al (2005) but without any heuristics
		- "event-pair": 					event-pair distance (e.g. as used in Schulmeister and Wheeler [2004])
}
  \item{seq_len}{
		The number of elements in the sequences (total number of elements in the data set, even if some are missing in these specific sequences)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	Variable, see description of 'output' argument above.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
	PGi uses a penalty system for accounting for unscorable data. Essentially this system assigns a penalty for the evolution "creation" or "loss" of an element. While biologically, this doesn't necesarily make sense, it is required in this implementation. Final tree lengths as reported by PGi correct for this and only return the number of sequence heterochronies.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{

	anc.seq<-matrix(c(1,2,2,2,2,3),1,6)
	der.seq<-matrix(c(1,2,2,2,3,3),1,6)
	print(edit.cost(anc.seq,der.seq,output=2,method="parsimov",seq_len=6))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
