\name{metabolites}
\alias{metabolites}
\title{Identifies the list of unique metabolites for a set of stoichiometric reactions
}
\description{
This function identifies the list of unique metabolites (reactants and products) for a set of stoichiometric reactions. This function assumes that all your reactions have the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
\usage{
metabolites(reactionList,woCompartment)
}
\arguments{
  \item{reactionList}{ A set of stoichiometric reactions with the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]}
  \item{woCompartment}{ A logical value (TRUE/FALSE) if compartments should be removed}
}
\examples{
# Loading the stoichiometric reactions from the glutamine/glutamate cycle
data(glugln)

# Extracting the metabolites from a set of reactions
metabolites(glugln)
# ...
# [272] "Quinone[x]"
# [273] "Quinone[c_n]"
# [274] "Quinone[b]"
# [275] "Succinate[x]"
# [276] "Succinate[b]"

# Extracting metabolites from a unique reaction
reaction <- glugln[2]
metabolites(reaction)
# [1] "alpha-D-Glucose 6-phosphate[r_n]" "H2O[r_n]" "alpha-D-Glucose[r_n]" "Orthophosphate[r_n]"

# Extracting the metabolites from a set of reactions without compartments
metabolites(glugln, woCompartment = TRUE)
# ...
# [64] "hydrogencarbonate"                                                              
# [65] "(S)-Lactate"                                                                    
# [66] "D-glutamic acid"                                                                
# [67] "FAD"                                                                            
# [68] "Ferricytochrome c2"                                                             
# [69] "Hydroquinone" 
}
\keyword{ metabolites }
\keyword{ stoichiometric }
