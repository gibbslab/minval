#' @export products
#' @title Identify the products of a stoichometric reaction
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @description This function identifies the products for a set of stoichometric reactions.
#' @param reactionList A set of stoichiometric reaction with the following format: 
#' 
#' \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"} 
#' 
#' Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}". 
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#' @return A vector with the identified products in the reaction, or a list if a set of stoichiometric reactions was given.
#' @examples
#' #' # Loading data
#' glycolysisKEGG <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
#'
#' # Extracting products for a single reaction
#' products(reactionList = "ADP[c] + Phosphoenolpyruvate[c] => Pyruvate[c] + ATP[c]")
#' 
#' # Extracting products for a set of stoichiometric reactions
#' products(reactionList = glycolysisKEGG$REACTION)
#' 
products <- function(reactionList){
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Remove reaction with invalid syntax
  reactionList <- reactionList[validateSyntax(reactionList)]
  # Extract reactants for irreversible reactions
  reaction <- strsplit(reactionList,"[[:blank:]]+=>[[:blank:]]+")
  reaction[lengths(reaction)>1] <- lapply(reaction[lengths(reaction)>1],function(reaction){reaction[[2]]})
  # Extract metabolites for reversible reactions
  reaction <- strsplit(unlist(reaction), "[[:blank:]]+<=>[[:blank:]]+")
  # Split independient reactants
  reaction <- lapply(reaction, function(reaction){strsplit(reaction,"[[:blank:]]+\\+[[:blank:]]+")})
  # Remove spaces and report uniques
  reaction <- lapply(reaction, function(reaction){unique(removeSpaces(unlist(reaction)))})
  # Use a regex to extract stoichiometric coefficients and separate the metabolite name
  products <- lapply(reaction, function(reaction){unique(removeCoefficients(reaction))})
  if (length(products)==1){
    return(unlist(products))
  } else {
    return(products)
  }
}

