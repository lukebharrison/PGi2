\name{pgi.to.phylo}
\alias{pgi.to.phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO CONVERT A PGI DATA OBJECT TO AN APE::PHYLO OBJECT
}
\description{
	The PGi data object is used to store the phylogenetic topology and sequences as well as execution parameters and results. To faciliate plotting and export, this function converts the PGi data object into a "phylo" object that can be used with ape. This funciton add several fields to the vanilla "phylo" object including, for example, the sequence heterochronies.
}
\usage{
pgi.to.phylo(pgi.tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pgi.tree}{
	A PGi data object of class "pgi.tree".
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
	vel.phylo<-pgi.to.phylo(velhagen)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
