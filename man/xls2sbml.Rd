\name{xls2sbml}
\alias{xls2sbml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Converts XLS reconstruction to a SBML file
%%  ~~function to do ... ~~
}
\description{
A function to convert XLS reconstruction to a SBML file}
\usage{
xls2sbml(infile,outfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{A XLS file with the following columns: \code{'ID'}, \code{'EQUATION'}, \code{'LOWER.BOUND'}, \code{'UPPER.BOUND'} and \code{'OBJECTIVE'}}
  \item{outfile}{A writeable path for a XML file}
}
\note{
This function includes a modified copy of saveSBML function included in SBMLR package. Original 'saveSBML' function was writed by Tomas Radivoyevitch. Please see: https://www.bioconductor.org/packages/release/bioc/html/SBMLR.html
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ xls }% use one of  RShowDoc("KEYWORDS")
\keyword{ sbml }% __ONLY ONE__ keyword per line
