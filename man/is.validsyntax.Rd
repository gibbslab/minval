\name{is.validsyntax}
\alias{is.validsyntax}
\title{ Evaluates if a stoichiometric reaction has a valid syntax
}
\description{
This function evaluates if a stoichiometric reaction has a valid syntax. It assumes that reaction have the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
\usage{
is.validsyntax(reaction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reaction}{ A stoichiometric reaction with the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
}
\examples{
## For a specific reaction
is.validsyntax("ATP[m_n] + beta-D-Glucose[m_n] => ADP[m_n] + beta-D-Glucose 6-phosphate[m_n]")
# [1] TRUE

## For a set of metabolites
# Loading a set of stoichiometric reactions from the glutamate/glutamine cycle
data("glugln")

# Validating syntaxis
sapply(glugln[1:10],is.validsyntax)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
