\name{chebi.id}
\alias{chebi.id}
\title{ Returns the ChEBI id asociated to a compound name
}
\description{
This function takes a compound name as an input and returns their associated ChEBI id
}
\usage{
chebi.id(metabolite)
}
\arguments{
  \item{metabolite}{A compound name
}
}
\examples{
## For a specific metabolite
chebi.id("ATP")

## For a set of metabolites
mets<-metabolites(glugln[1:3],woCompartment = TRUE)
sapply(mets, chebi.id, simplify = FALSE)
# $`alpha-D-glucose 1-phosphate`
# [1] "29042"

# $`alpha-D-Glucose 6-phosphate`
# [1] "17665"

# $water
# [1] "15377"

# $ATP
# [1] "15422"

# $`alpha-D-Glucose`
# [1] "17925"

# $ADP
# [1] "16761"

# $`phosphate(3-)`
# [1] "18367"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
