\name{summary.multi.pgi.tree}
\alias{summary.multi.pgi.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A FUNCTION TO SUMMARIZE THE CONTENT OF A MULTI PGI DATA OBJECT
}
\description{
summary.multi.pgi.tree() summarizes the content of a multi PGi data object (a list of PGi trees). This function summaries information on type (e.g. empty, filled, consensus) and tree lengths, if available.
}
\usage{
summary.multi.pgi.tree(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	An object of class "multi.pgi.tree", which is a simple list of "pgi.tree" objects.
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
	\code{\link{plot.pgi.tree}},\code{\link{summary.pgi.tree}}
}
\examples{

	data(velhagen)
	data(sanchez.villagra2002) 
	multi.t<-list(velhagen,sanchez.villagra2002)
	class(multi.t)<-"multi.pgi.tree"
	summary.multi.pgi.tree(multi.t)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
