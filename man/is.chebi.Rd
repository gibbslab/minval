\name{is.chebi}
\alias{is.chebi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Evaluates if a compound name is a ChEBI name
}
\description{
This function takes a compound name and returns a logical value (TRUE) if it is a ChEBI name}
\usage{
is.chebi(metabolite)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{metabolite}{A compound name
}
}
\examples{
## For a specific metabolite
is.chebi("ATP")
# [1] TRUE

## For a set of metabolites
mets<-metabolites(glugln[1:3],woCompartment = TRUE)
sapply(mets,is.chebi,simplify = FALSE)
# $`alpha-D-glucose 1-phosphate`
# [1] TRUE

# $`alpha-D-Glucose 6-phosphate`
# [1] TRUE

# $water
# [1] TRUE

# $ATP
# [1] TRUE

# $`alpha-D-Glucose`
# [1] TRUE

# $ADP
# [1] TRUE

# $`phosphate(3-)`
# [1] TRUE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
