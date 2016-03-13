\name{orphan.reactants}
\alias{orphan.reactants}
\title{Identifies the orphan reactants for a set of stoichometric reactions
}
\description{
This function identifies the orphan reactants (reactants never produced in any reaction) for a set of stoichiometric reactions. Possible candidates to be introduced into the system by exchange reactions or by adding more internal reactions. This function assumes that all your reactions have the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
\usage{
orphan.reactants(reactionList,byCompartment)
}
\arguments{
  \item{reactionList}{ A set of stoichiometric reactions}
  \item{byCompartment}{A logical value (TRUE/FALSE) if orphans should be reported by compartment}
}
\examples{
# Loading a set of stoichiometric reactions from the glutamate/glutamine cycle
data("glugln")

# Identifying orphan reactants from a sample of 100 stoichiometric reactions (*)
orphan.reactants(glugln)
# [1] "alpha-D-Glucose 6-phosphate[r_n]"    
# [2] "water[r_n]"                          
# [3] "2,3-bisphospho-D-glyceric acid[r_n]"
# [4] "GTP[c_n]"                            
# [5] "oxaloacetic acid[m_n]"
# ...
# [28] "diphosphate(4-)[m_a]"
}
\keyword{ orphan }
\keyword{ reactants }
\keyword{ metabolites }
