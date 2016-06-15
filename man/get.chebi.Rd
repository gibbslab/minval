\name{get.chebi}
\alias{get.chebi}
\title{Returns associated data for a compound name from the ChEBI database
}
\description{This function takes a compound or a set of compounds names as an input and returns associated data from the ChEBI database.
}
\usage{
get.chebi(metabolite,get)
}
\arguments{
  \item{metabolite}{A compound or a set of compounds names}
  \item{get}{A character string specifying the associated data to be retrieved, must be one of \code{"id"} }
}
\examples{
## For a specific metabolite
get.chebi("ATP")
}
\keyword{ ChEBI }% use one of  RShowDoc("KEYWORDS")
\keyword{ data }% __ONLY ONE__ keyword per line
