\name{compartments}
\alias{compartments}
\title{Identifies the compartments for a set of metabolites
}
\description{
This function identifies the compartments for a set of compartmentalized metabolites. This function assumes that all your reactions have the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
\usage{
compartments(metabolites)
}
\arguments{
  \item{metabolites}{ A set of compartimentalized metabolites
}
}
\examples{
# Loading a set of stoichiometric reactions from the glutamate/glutamine cycle
data("glugln")

# For a set of metabolites from stoichiometric reactions
compartments(metabolites(glugln))

# For a specific metabolite
compartments("NADH[c]")
}
\keyword{ metabolite }
\keyword{ compartments }
