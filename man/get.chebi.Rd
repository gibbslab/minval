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
  \item{get}{A character string specifying the associated data to be retrieved, must be one of \code{"id"}, \code{"formula"} or {"candidates"} }
}
\examples{
## Retrieve ChEBI id for a specific metabolite name
get.chebi("ATP",get = "id")
# "15422" 

## Retrieve molecular formula for a specific metabolite name
get.chebi("ATP",get = "formula")
# "C10H16N5O13P3"

## Retrieve candidates for a specific metabolite name
get.chebi("ATP",get = "candidates")
# [1] "atp"                                                                                           
# [2] "datp"                                                                                          
# [3] "5-phosphoribosyl-atp"                                                                          
# [4] "atp synthase inhibitor"                                                                        
# [5] "atp-sugar"               
# .
# .
# .
}
\keyword{ ChEBI }
\keyword{ data }
