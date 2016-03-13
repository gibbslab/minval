\name{products}
\alias{products}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identifies the products for a stoichometric reaction
}
\description{
This function identifies the products for a stoichometric reaction. This function assumes that all your reactions have the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
\usage{
products(reaction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reaction}{ A stoichiometric reaction with the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
%%     ~~Describe \code{rxn} here~~
}
}
\examples{
# For a stoichiometric reaction
reaction <- "ADP[c_n] + Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n]"

# Extract the products
products(reaction)
# [1] "ATP[c_n]"      
# [2] "Pyruvate[c_n]"

# If stoichiometric reaction is reversible, function returns all metabolites
products("ADP[c_a] + ATP[m_a] <=> ADP[m_a] + ATP[c_a]")
# [1] "ADP[c_a]"
# [2] "ATP[m_a]"
# [3] "ADP[m_a]"
# [4] "ATP[c_a]"
}
\keyword{ products }
\keyword{ metabolites }% __ONLY ONE__ keyword per line
