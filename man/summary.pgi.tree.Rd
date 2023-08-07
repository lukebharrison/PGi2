\name{summary.pgi.tree}
\alias{summary.pgi.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A FUNCTION TO SUMMARIZE THE CONTENT OF A PGI DATA OBJECT
}
\description{
summary.pgi.tree() summarizes the content of a PGi data object. Basic parameters including sequence length, number of nodes and the phylogenetic topology are always printed. By default, only the first two developmental sequence are printed, this behaviour is overriden using print.all.seqs=TRUE. In the case of filled and consensus trees, tree length and consensus information is also printed. If PGi execution parameters have been stored in the data object, they are also printed. 
}
\usage{
summary.pgi.tree(object, print.all.seqs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	An onject of class "pgi.tree" of any type (empty, filled, consensus, or superconsensus).
}
  \item{print.all.seqs}{
	If (TRUE), then print all developmental sequences, if (FALSE, default) print only the first two sequences.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	None.
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
	\code{\link{plot.pgi.tree}},\code{\link{summary.multi.pgi.tree}}
}
\examples{

	data(velhagen)
	summary.pgi.tree(velhagen)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
