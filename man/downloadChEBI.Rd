% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downloadChEBI.R
\name{downloadChEBI}
\alias{downloadChEBI}
\title{Download the ChEBI database}
\usage{
downloadChEBI(release = "latest", woAssociations = FALSE)
}
\arguments{
\item{release}{A character string with the release number of the ChEBI database version to be downloaded, by default the \code{'latest'} release is downloaded.}

\item{woAssociations}{A logical value \code{'TRUE'} or \code{'FALSE'} if a light version of the ChEBI database without associations should be returned.}
}
\value{
A data.frame with the following data associated to the ChEBI compounds: \itemize{
\item \code{'ID'}: The unique identifer
\item \code{'ChEBI'}: The name recommended for use in biological databases
\item \code{'KEGG'}: The associated name(s) in the KEGG database
\item \code{'IUPAC'}: The name(s) generated according to recommendations of IUPAC
\item \code{'MetaCyc'}: The associated name(s) in the MetaCyc database
\item \code{'ChEMBL'}: The associated name(s) in the ChEMBL database
\item \code{'FORMULA'}: The molecular formula
\item \code{'MASS'}: The molecular mass
\item \code{'MONOISOTOPIC'}: The molecular monoisotopic mass
\item \code{'CHARGE'}: The molecular net charge
}
If woAssociations is \code{'TRUE'} a A data.frame with the following data is returned: \itemize{
\item \code{'NAME'}: The name(s) associated in several biological databases
\item \code{'FORMULA'}: The molecular formula
\item \code{'MASS'}: The molecular mass
\item \code{'MONOISOTOPIC'}: The molecular monoisotopic mass
\item \code{'CHARGE'}: The molecular net charge
}
}
\description{
This function downloads the compounds, formulas, masses and charges from the selected release of the ChEBI database.
The ChEBI database (Chemical Entities of Biological Interest), is a database and ontology of molecular entities focused on 'small' chemical compounds.
}
\examples{
\dontrun{
# Download ChEBI database with associations
ChEBI <- downloadChEBI(release = '142')

# Download ChEBI database without associations
ChEBI <- downloadChEBI(release = '142', woAssociations = TRUE)
 }
}
\seealso{
The ChEBI database webpage: https://www.ebi.ac.uk/chebi/
}
\author{
Daniel Camilo Osorio <dcosorioh@tamu.edu>
}
