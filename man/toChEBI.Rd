\name{toChEBI}
\alias{toChEBI}
\title{ Translates compounds names to ChEBI ids in a stoichiometric reaction
}
\description{ This function translates compounds names to ChEBI ids or molecular formulas in a stoichiometric reaction}
\usage{
toChEBI(reaction, formula = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reaction}{A stoichiometric reaction}
  \item{formula}{ A logical value, if TRUE translates the compounds names to molecular formulas}
}
\examples{
## For a specific reaction
toChEBI("ADP[c_n] + Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n]")
# [1] "1 16761 + 1 18021 => 1 15422 + 1 15361"

toChEBI("ADP[c_n] + Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n]",formula = TRUE)
# [1] "1 C10H15N5O10P2 + 1 C3H4O6P => 1 C10H16N5O13P3 + 1 C3H3O3"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
