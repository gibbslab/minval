#' @export metabolites
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Identify the list of metabolites for a set of stoichiometric reactions
#' @description This function identifies the list of metabolites for a set of stoichiometric reactions. If \code{'woCompartment'} is \code{'TRUE'} compartment label is removed. If \code{'uniques'} is \code{'TRUE'}, list of uniques is returned.
#' @param reactionList A set of stoichiometric reaction with the following format: 
#' 
#' \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"} 
#' 
#' Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}". 
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#' @param woCompartment A boolean value \code{'TRUE'} or \code{'FALSE'} to indicate if compartment label should be removed
#' @param uniques A boolean value \code{'TRUE'} or \code{'FALSE'} to indicate if uniques must be returned
#' @return A list of metabolites for a set of stoichiometric reactions

metabolites <- function(reactionList, woCompartment = FALSE, uniques=TRUE){
  # Split reactions by arrow symbol
  reaction <- strsplit(as.vector(reactionList),"[[:blank:]]+<?=>[[:blank:]]*")
  # Split sections by plus symbol
  reaction <- lapply(reaction, function(reaction){strsplit(unlist(reaction),"[[:blank:]]+\\+[[:blank:]]+")})
  # Remove extra spaces around metabolite names
  reaction <- lapply(reaction, function(reaction){removeSpaces(unlist(reaction))})
  # Remove coefficients
  reaction <- lapply(reaction, function(reaction){removeCoefficients(reaction)})
  # Convert metabolite list to a vector
  metabolites <- unlist(reaction)
  # Option to remove compartments
  if (woCompartment == TRUE) {
    metabolites <- removeCompartment(metabolites)
  }
  # Option to return uniques
  if (uniques == TRUE) {
    metabolites <- unique(metabolites)
  }
  # Return a metabolite list
  return(metabolites)
}
