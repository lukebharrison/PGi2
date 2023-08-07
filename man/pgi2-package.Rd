\name{pgi2-package}
\alias{pgi2-package}
\alias{pgi2}
\docType{package}
\title{
PGi2: Parsimov-scored Genetic Inference Method for the Analysis of Developmental Sequences and Sequence Heterochrony
PGI2
}
\description{
This package implements the PGi algorithm described in Harrison and Larsson (2008), Syst. Biol. 57:378-387. This package is used in the evolutionary analysis of developmental sequences on phylogenies. Please see the README included with this package, and the vignettes for a demonstration of its use. Please also find documentaiton for each individual function in the R documentation system. The main function is called pgi() and ?pgi is a good place to start for documentation.
}
\details{
\tabular{ll}{
Package: \tab pgi2\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-02-07\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
Please refer to the vignettes and the Harrison and Larsson (2008) article for a full description of usage of the PGi package.
}
\author{
Luke B. Harrison

Maintainer: Luke Harrison <luke.harrison@mail.mcgill.ca>
}
\references{
Harrison, Luke B. and Larsson, Hans C. E. (2008) 'Estimating Evolution of Temporal Sequence Changes: A Practical Approach to Inferring Ancestral Developmental Sequences and Sequence Heterochrony', Systematic Biology, 57:3, 378—387.
}
\keyword{
	package
	heterochrony
	phylogeny
}
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
	# Fully interactive execution
	data(velhagen)
	pgi(velhagen,interactive=TRUE)

	# Non-interactive execution
	data(sanchez.villagra.2002)
	sv2002.con.trees<-pgi(sanchez.villagra.2002,nruns=c(2,2),inf.params=list(heuristic="pgi",cycles=100,replicates=100,ret.anc.seq=100,edit.cost.func="parsimov",simultaneity=TRUE),con.params=list(con.type="semi.exhaustive",semi.ex.con.max.n=5000,edit.cost.func="parsimov"),verbosity=1)
	sv2002.supercon<-pgi.supercon(sv2002.con.trees,supercon.tol=0)
}

