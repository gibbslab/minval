\name{reactants}
\alias{reactants}
\title{Identifies the reactants for a stoichometric reaction
}
\description{
This function identifies the reactants for a stoichometric reaction. This function assumes that all your reactions have the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
\usage{
reactants(rxn)
}
\arguments{
  \item{rxn}{A stoichiometric reaction with the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
}
\examples{
# For a stoichiometric reaction
reaction <- "ADP[c_n] + Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n]"

# Extract the products
reactants(reaction)
# [1] "ADP[c_n]"                 
# [2] "Phosphoenolpyruvate[c_n]"

# If stoichiometric reaction is reversible, function returns all metabolites
reactants("ADP[c_a] + ATP[m_a] <=> ADP[m_a] + ATP[c_a]")
# [1] "ADP[c_a]"
# [2] "ATP[m_a]"
# [3] "ADP[m_a]"
# [4] "ATP[c_a]"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ reactants }
\keyword{ metabolites }% __ONLY ONE__ keyword per line
