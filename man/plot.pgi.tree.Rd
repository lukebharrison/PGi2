\name{plot.pgi.tree}
\alias{plot.pgi.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	A FUNCTION TO PLOT A PGI OBJECT
}
\description{
	plot.pgi.trees() plots a PGi data object. In the case of empty PGi objects, the topology and optionally, the sequences at the tips are plotted. In the case of consensus objects, the topology, branch lengths, sequence heterochronies, ancestral sequences and heterochrony support values can all be plotted (see argument details below). Note: this function relies on the ape package for plotting functions. In the case of very large data sets, it is important to select a large enough plot area or 
}
\usage{
plot.pgi.tree(pgi.tree, seq.het.thres = 0.5, print.support = FALSE, show.anc.seq = FALSE, show.tip.seq = FALSE, show.bls = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pgi.tree}{
	PGi data object (of class "pgi.tree" or "multi.pgi.tree") of any type
}
  \item{seq.het.thres}{
	What pseudoconensus support threshold to use when printing recovered heterochronies (e.g. 0.5 for majority-rule, 1 for strict).
}
  \item{print.support}{
	TRUE/FALSE: Should pseudoconsensus support values be printed parenthetically after recovered sequence heterochronies.
}
  \item{show.anc.seq}{
	TRUE/FALSE: Should pseudoconsensus ancestral sequences be printed at the nodes of the phylogeny (can be messy with large data sets / long sequence lengths).
}
  \item{show.tip.seq}{
	TRUE/FALSE: Should the developmental seqeunces at the tips be printed on teh phylogeny (Can be messy with long sequences).
}
  \item{show.bls}{
	TRUE/FALSE: should branches be drawn proportionally to the number of sequence heterochronies on each branch? Note branch lenghts correspond to the number of heterochronies under the threshold inputted in seq.het.thres or 0.5 by default.
}
  \item{\dots}{
	Additional arguments to either ape::plot.phylo (on which this function is build) or to plot().
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
		

	data(velhagen)
	plot.pgi.tree(velhagen, show.tip.seq=TRUE)
	data(velhagen.pseudocon)
	plot(velhagen.pseudocon,show.anc.seq=TRUE,seq.het.thres=0.5,print.support=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
