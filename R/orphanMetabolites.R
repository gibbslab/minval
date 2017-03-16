#' @export orphanMetabolites
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Identify the orphan metabolites of a set of stoichiometric reactions
#' @description This function identifies the orphan metabolites (metabolites not produced or not consumed in any other reaction or just involved in one reaction) for a set of stoichometric reactions.
#' @param reactionList A set of stoichiometric reaction with the following format:
#'
#' \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"}
#'
#' Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}".
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#' @param actingAs A text string that specifies the type of metabolite to be returned; only \code{'reactant'} and \code{'product'} are supported.
#' @param byCompartment A boolean value \code{'TRUE'} or \code{'FALSE'} to indicate if orphan reactants should be reported by compartment
#' @return If \code{byCompartment == FALSE}, a vector with orphan reactants is returned, in opposite case a list is returned.
#' 
#' If \code{actingAs == 'reactant'}, metabolites not produced in any other reaction or just are involved in one reaction are returned.
#' 
#' If \code{actingAs == 'products'}, metabolites not consumed in any other reaction or just are involved in one reaction are returned.
#' @examples
#' # Loading data
#' glycolysisKEGG <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"), sep = "\t")
#' 
#' # Extracting orphan reactants
#' orphanMetabolites(reactionList = glycolysis$REACTION, actingAs = "reactant")
#'
#' # Extracting orphan products by compartment
#' orphanMetabolites(reactionList = glycolysis$REACTION, actingAs = "product", byCompartment = TRUE)


orphanMetabolites <- function(reactionList, actingAs = NULL, byCompartment = FALSE){
  # Convert to a vector
  reactionList <- as.vector(x = reactionList)
  # Remove reaction with invalid syntax
  reactionList <- reactionList[validateSyntax(reactionList = reactionList)]
  # Extract all reactants and products
  reactant <- unique(unlist(reactants(reactionList = reactionList)))
  product <- unique(unlist(products(reactionList = reactionList)))
  # Identify metabolites with frequency equal to 1
  met <- table(metabolites(reactionList = reactionList, uniques = FALSE))
  met <- names(met)[met==1]
  # Identify orphans
  if (is.null(actingAs)){
    orphan <- unique(c(reactant[!reactant %in% product], product[!product %in% reactant], met))
  } else if (actingAs == "reactant"){
    orphan <- unique(c(reactant[!reactant %in% product], met[met %in% reactant]))
  } else if (actingAs == "product"){
    orphan <- unique(c(product[!product %in% reactant], met[met %in% product]))
  } else {
    stop("actingAs should be one of: 'reactant' or 'product'")
  }
  # Return
  if (byCompartment == TRUE) {
    orphan <- sapply(compartments(orphan), function(compartment){orphan[grep(paste0("\\[", compartment, "\\]"), orphan)]}, simplify = FALSE)
  }
  return(orphan)
}
