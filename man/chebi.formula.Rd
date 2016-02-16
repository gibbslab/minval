\name{chebi.formula}
\alias{chebi.formula}
\title{ Returns the molecular formula associated to a ChEBI compound name
}
\description{
This function takes a compound or a set of compounds names as an input and returns a list with molecular formulas.
}
\usage{
chebi.formula(metabolite)
}
\arguments{
  \item{metabolite}{A compound or a set of compounds names
}
}
\examples{
## For a specific metabolite
chebi.formula("ATP")
# [1] "C10H16N5O13P3"

## For a set of metabolites
mets<-metabolites(glugln[1:3],woCompartment = TRUE)
sapply(mets,chebi.formula,simplify = FALSE)

# $`alpha-D-glucose 1-phosphate`
# [1] "C6H13O9P"

# $`alpha-D-Glucose 6-phosphate`
# [1] "C6H13O9P"

# $water
# [1] "H2O"

# $ATP
# [1] "C10H16N5O13P3"

# $`alpha-D-Glucose`
# [1] "C6H12O6"

# $ADP
# [1] "C10H15N5O10P2"

# $`phosphate(3-)`
# [1] "O4P"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
