#' @export orphanProducts
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Identify the orphan products of a set of stoichometric reactions
#' @description This function identifies the orphan products (products that are not consumed in any other reaction or just are involved in one reaction) for a set of stoichometric reactions.
#' @param reactionList A set of stoichiometric reaction with the following format: \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"} Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}". 
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#' @param byCompartment A boolean value \code{'TRUE'} or \code{'FALSE'} to indicate if orphan products must be reported by compartment
#' @return If \code{byCompartment == FALSE}, a vector with orphan products is returned, in opposite case a list is returned.
#' @examples 
#' # Loading data
#' glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))

#' # Removing stoichiometric reactions without valid syntax
#' glycolysis <- mapReactions(
#' reactionList = isValidSyntax(glycolysis$REACTION),
#' referenceData = glycolysis,
#' by = "bool"
#' )
#' 
#' # Extracting orphan products
#' orphanProducts(reactionList = glycolysis$REACTION, byCompartment = FALSE)
#' 
#' # Extracting orphan by compartment
#' orphanProducts(reactionList = glycolysis$REACTION, byCompartment = TRUE)
orphanProducts <- function(reactionList, byCompartment=FALSE){
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Remove reaction with invalid syntax
  reactionList <- reactionList[isValidSyntax(reactionList)]
  # Extract all reactants
  reactant <- unique(unlist(reactants(reactionList)))
  # Extract all products
  product <- unique(unlist(products(reactionList)))
  # Possible candidates to be introduced into the system by exchange reactions or by adding more internal reactions.
  orphan <- rowSums(stoichiometricMatrix(reactionList)!=0)
  orphan <- c(names(orphan)[orphan<2],product[!(product%in%reactant)])
  orphan <- unique(orphan[!orphan%in%reactant[!reactant%in%product]])
  if(length(orphan)==0){
    return (NA)
  } else {
    if (byCompartment == TRUE){
      # Return orphans by compartment
      sapply(compartments(orphan), function(compartment){orphan[grep(paste0("\\[",compartment,"\\]"),orphan)]}, simplify = FALSE)
    } else {
      # Return all reactants never produced in any reaction.
      return(orphan)
    }
  }
}



