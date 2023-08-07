\name{extract.seqs.from.tree}
\alias{extract.seqs.from.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	UTILITY FUNCTION TO EXTRACT ALL DEVELOPMENTAL SEQUENCES FROM A PGI DATA OBJECT (TREE)
}
\description{
	This function is used on a pgi data object's tree field (E.g. pgi.tree$tree) and is used to extract all tip developmental sequences.
}
\usage{
extract.seqs.from.tree(tree, seq_len, seqs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
		A raw PGi tree data structure (found in pgi.tree$tree)
}
  \item{seq_len}{
		The sequence length (number of elements), found in pgi.tree$seq.len
}
  \item{seqs}{
		Names of specific sequences to recover, if possible. Defaults to all sequences.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A matrix of ranked developmental sequences (rows) and tipnames as rownames.
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
extracted.sequences<-extract.seqs.from.tree(velhagen$tree,velhagen$seq.len)	


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
