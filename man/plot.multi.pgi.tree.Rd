\name{plot.multi.pgi.tree}
\alias{plot.multi.pgi.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	A FUNCTION TO PLOT A MULTI PGI OBJECT
}
\description{
	plot.multi.pgi.trees() plots a set of PGi data objects, by calling plot.pgi.tree(). The format of the plot is set up using the format arugment, trees can be plotted in a matrix in a single plot or in multiple individual plots.
}
\usage{
plot.pgi.tree(pgi.tree, format=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multi.pgi.tree}{
	PGi data objects (of class "multi.pgi.tree") of any type
}
  \item{format}{
	Plot format: set the format of the plot to be generated of the multiple pgi.tree objects. This will be passed to par(mfrow=c(X,Y)), and should be specified as a two integer vector. If no format is provided, one will be calculated based on the number of trees. If the numbers of the trees is high, this will lead to plots that are too small.
}
  \item{...}{
	Arguments to be passed to plot.pgi.tree() for plotting options for the individual plots. Please see the plot.pgi.tree() man page for details.
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
	velhagen.con.trees<-pgi(velhagen,nruns=c(8,1),inf.params=list(heuristic="pgi",cycles=100,replicates=100,ret.anc.seq=100,edit.cost.func="parsimov",simultaneity=TRUE),con.params=list(con.type="semi-exhaustive",semi.ex.con.max.n=5000,edit.cost.func="parsimov"),verbosity=1)
	plot.multi.pgi.tree(velhagen.con.trees)
	

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
