\name{stoichiometric.matrix}
\alias{stoichiometric.matrix}
\title{Returns a stoichiometry matrix for a set of stoichiometric reactions
}
\description{ This function returns a stoichiometry matrix where each row corresponds to a metabolite, and each column corresponds to a reaction.
}
\usage{
stoichiometric.matrix(reactionList)
}
\arguments{
  \item{reactionList}{ A set of stoichiometric reactions with the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]}
}

\examples{
# loading data
data(glugln)

# create stoichiometric matrix
stoichiometric.matrix(glugln)
}
\keyword{ stoichiometric }
\keyword{ matrix }
