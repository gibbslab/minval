\name{metabolites}
\alias{metabolites}
\title{Identifies the list of unique metabolites for a set of stoichiometric reactions
}
\description{
This function identifies the list of unique metabolites (reactants and products) for a set of stoichiometric reactions. This function assumes that all your reactions have the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
\usage{
metabolites(reactionList,woCompartment,uniques)
}
\arguments{
  \item{reactionList}{ A set of stoichiometric reactions with the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]}
  \item{woCompartment}{ A logical value (TRUE/FALSE) if compartments should be removed}
  \item{uniques}{ A logical value (TRUE/FALSE) if unique metabolites should be reported}
}
\examples{
# Loading the stoichiometric reactions from the glutamine/glutamate cycle
data(glugln)

# Extracting the metabolites from a set of reactions
metabolites(glugln)
# ...
# [201] "oxygen atom[m_a]"                                                                    
# [202] "Ferrocytochrome c2[m_a]"                                                             
# [203] "Ferricytochrome c2[m_a]"                                                             
# [204] "diphosphate(4-)[m_a]"                                                                
# [205] "(R)-Lactate[m_a]"  

# Extracting metabolites from a unique reaction
reaction <- glugln[2]
metabolites(reaction)
# [1] "alpha-D-Glucose 6-phosphate[r_n]" "H2O[r_n]" "alpha-D-Glucose[r_n]" "phosphate(3-)[r_n]"

# Extracting the metabolites from a set of reactions without compartments
metabolites(glugln, woCompartment = TRUE)
# ...
# [64] "oxygen atom"                                                                    
# [65] "Ferrocytochrome c2"                                                             
# [66] "Ferricytochrome c2"                                                             
# [67] "diphosphate(4-)"                                                                
# [68] "(R)-Lactate"                                                                    
# [69] "hydrogencarbonate"    
}
\keyword{ metabolites }
\keyword{ stoichiometric }
