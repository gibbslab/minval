\name{orphan.products}
\alias{orphan.products}
\title{Identifies the orphan products for a set of stoichometric reactions
}
\description{
This function identifies the orphan products (products that are not consumed in any other reaction. Possible dead-ends) for a set of stoichometric reactions. This function assumes that all your reactions have the following format L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
\usage{
orphan.products(reactionList,byCompartment)
}
\arguments{
  \item{reactionList}{A set of stoichiometric reactions}
  \item{byCompartment}{A logical value (TRUE/FALSE) if orphans should be reported by compartment}
}

\examples{
# Loading a set of stoichiometric reactions from the glutamate/glutamine cycle
data("glugln")

# Identifying dead-ends from a sample of 100 stoichiometric reactions
orphan.products(glugln)
# [1] "alpha-D-Glucose[r_n]"           
# [2] "phosphate(3-)[r_n]"             
# [3] "2-phospho-D-glyceric acid[r_n]"
# [4] "GDP[c_n]"                       
# ...
# [20] Ferricytochrome c2[m_a]"

}
\keyword{ orphan }
\keyword{ products }
\keyword{ dead-ends }
