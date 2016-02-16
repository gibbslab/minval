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
unbalanced("ADP[c_n] + Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n] + Water[c_n]")
# [1] TRUE

unbalanced("Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n]",show.formulas = TRUE)
#      [,1]                                                                          
# [1,] "ADP[c_n] + Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n] + Water[c_n]"
#      [,2]                                                               
# [1,] "1 C10H15N5O10P2 + 1 C3H4O6P => 1 C10H16N5O13P3 + 1 C3H3O3 + 1 H2O"

## For a set of reactions
# Loading the stoichiometric reactions from the glutamine/glutamate cycle
data(glugln)

sapply(glugln[1:5],unbalanced,simplify = FALSE)
# $`alpha-D-glucose 1-phosphate[c_n] <=> alpha-D-Glucose 6-phosphate[c_n]`
# [1] FALSE

# $`alpha-D-Glucose 6-phosphate[r_n] + water[r_n] => alpha-D-Glucose[r_n] + phosphate(3-)[r_n]`
# [1] TRUE

# $`ATP[m_n] + alpha-D-Glucose[m_n] <=> ADP[m_n] + alpha-D-Glucose 6-phosphate[m_n]`
# [1] FALSE

# $`ATP[m_n] + beta-D-Glucose[m_n] => ADP[m_n] + beta-D-Glucose 6-phosphate[m_n]`
# [1] FALSE

# $`alpha-D-Glucose[m_n] + ADP[m_n] <=> alpha-D-Glucose 6-phosphate[m_n] + AMP[m_n]`
# [1] FALSE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
