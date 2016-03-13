\name{unbalanced}
\alias{unbalanced}
\title{Evaluates if a stoichiometric reaction is mass-balanced
}
\description{
This function evaluates if a stoichiometric reaction is mass-unbalanced. if TRUE, an option to return molecular formula  is included}
\usage{
unbalanced(reaction,show.formulas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reaction}{A stoichiometric reaction with the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]}
  \item{show.formulas}{A logical value (TRUE/FALSE) if molecular formula should be showed for unbalanced reactions}
}
\examples{
## For a specific reaction
unbalanced("Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n] + Water[c_n]")
# [1] TRUE

unbalanced("Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n]", show.formulas = TRUE)
#      [,1]                                                                          
# [1,] "Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n]"
#      [,2]                                                               
# [1,] "1 C3H4O6P => 1 C10H16N5O13P3 + 1 C3H3O3"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
