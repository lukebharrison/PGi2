\name{pgi.read.nexus}
\alias{pgi.read.nexus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	FUNCTION TO READ A NEXUS FORMATTED FILE INTO A PGI DATA OBJECT
}
\description{
	This function reads a NEXUS file which contains developmental sequences (see NEXUS file coding, in details, below). Note that this function requires the ape package and is based on the read.nexus() and read.nexus.data() functions from ape.
}
\usage{
pgi.read.nexus(file = NULL, interactive = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
		Filename of the NEXUS file containing developmental sequences and topology. Either fully qualified or simple filename in the current working directory.
}
  \item{interactive}{
		TRUE/FALSE: if TRUE and the filename is not provided, interactively asks for the file name or opens a TCLTK filechooser box to select the file.
}
}
\details{
NEXUS File Coding Details:
Enter the RANKED data as discrete characters. 
For example, if ABDC is the sequence of elements, enter: 
character 1 (A) = (rank) 1
character 2 (B) = (rank) 2
character 3 (C) = (rank) 4
character 4 (D) = (rank) 3 
and so on.

Enter unscorable (i.e. evolutionary loss of event) as 'Z'
Do not enter missing data (yet) - sorry: As a temporary fix, 
enter it as unscorable - biologically this doesnt really make
sense but it works out in the analysis.

For ranks between 10 and 35 (the max for now), use A,B,C,D, in the nexus file and so on.
It is possible to have 36+ but the the pgi tree must be manually edited - contact luke.harrison@mcgill.ca for details of how to do this.

Make sure the taxon names have no spaces and aren't too long (it'll copy those too). 
Enter a single phylogenetic topology and save it in to the nexus file in the standard parenthetical format (newick).

}
\value{
	A pgi.tree data object of $type = empty, with the $tree slot filled with the topology and sequence in the NEXUS file.
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
