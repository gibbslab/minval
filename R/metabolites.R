#' @export metabolites
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Identify the list of metabolites for a set of stoichiometric reactions
#' @description This function identifies the list of metabolites for a set of stoichiometric reactions. If \code{'woCompartment'} is \code{'TRUE'} compartment label is removed. If \code{'uniques'} is \code{'TRUE'}, list of uniques is returned.
#' @param reactionList 
#' @param woCompartment 
#' @param uniques 
#' @return A list of metabolites for a set of stoichiometric reactions
#' @examples 
#' # Loading data
#' glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
#' glycolysis <- mapReactions(reactionList = isValidSyntax(glycolysis$REACTION),referenceData = glycolysis,by = "bool")
#' 
#' # Extracting metabolites
#' metabolites(reactionList = glycolysis$REACTION)
#' 
#' # Extracting metabolites without compartments
#' metabolites(reactionList = glycolysis$REACTION, woCompartment = TRUE)
#' 
#' # Extracting redundant list of metabolites
#' metabolites(reactionList = glycolysis$REACTION, woCompartment = FALSE,uniques = FALSE)


metabolites <- function(reactionList, woCompartment = FALSE, uniques=TRUE){
  reaction <- strsplit(as.vector(reactionList),"[[:blank:]]+<?=>[[:blank:]]*")
  reaction <- lapply(reaction, function(reaction){strsplit(unlist(reaction),"[[:blank:]]+\\+[[:blank:]]+")})
  reaction <- lapply(reaction, function(reaction){.remove.spaces(unlist(reaction))})
  reaction <- lapply(reaction, function(reaction){.remove.coefficients(reaction)})
  metabolites <- unlist(reaction)
  if (woCompartment == TRUE) {
    metabolites <- .remove.compartment(metabolites)
  }
  if (uniques == TRUE) {
    metabolites <- unique(metabolites)
  }
  return(metabolites)
}
