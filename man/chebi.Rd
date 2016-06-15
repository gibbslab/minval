\name{chebi}
\alias{chebi}
\docType{data}
\title{ Public data from ChEBI database release 136
}
\description{
This dataset contains public data (61895 compounds) from the release 136 of the Chemical Entities of Biological Interest (ChEBI) database 
}
\usage{data("chebi")}
\format{
  A data frame with 61895 observations on the following 4 variables.
  \describe{
    \item{\code{id}}{the unique identifer}
    \item{\code{name}}{the name recommended for use in biological databases}
    \item{\code{candidates}}{other names together with an indication of their source}
    \item{\code{formula}}{Molecular formula}
  }
}
\source{
Dataset was compiled through the function \code{get.chebi.all} include in the RbioRXN package.
}
\references{
Hastings, J., de Matos, P., Dekker, A., Ennis, M., Harsha, B., Kale, N., Muthukrishnan, V., Owen, G., Turner, S., Williams, M., and Steinbeck, C. (2013) The ChEBI reference database and ontology for biologically relevant chemistry: enhancements for 2013. Nucleic Acids Res.
}
\examples{
data(chebi)
}
\keyword{datasets}
