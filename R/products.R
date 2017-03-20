#' @export products
#' @title Identify the products of a stoichometric reaction
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @description This function identifies the products for a set of stoichometric reactions.
#' @param reactionList A set of stoichiometric reaction with the following characteristics: \itemize{
#' \item Arrows symbols must be given in the form \code{'=>'} or \code{'<=>'}
#' \item Inverse arrow symbols \code{'<='} or other types as: \code{'-->'}, \code{'<==>'}, \code{'->'} will not be parsed and will lead to errors.
#' \item Arrow symbols and plus signs (\code{+}) must be surrounded by a space character
#' \item Stoichiometric coefficients must be surrounded by a space character and not by parentheses.
#' \item Each metabolite must have only one stoichiometric coefficient, substituents must be joined to metabolite name by a hyphen (\code{-}) symbol.
#' \item Exchange reactions have only one metabolite before arrow symbol
#' \item Compartments must be given between square brackets ([compartment]) joined at the end of metabolite name
#' }
#' Some examples of valid stoichiometric reactions are: \itemize{
#' \item \code{H2O[c] + Urea-1-Carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]}
#' \item \code{ADP[c] + Phosphoenolpyruvate[c] => ATP[c] + Pyruvate[c]}
#' \item \code{CO2[c] <=> }
#' }
#' @return A vector with the identified products in the reaction, or a list with the identified products in each reaction if a set of stoichiometric reactions was given.

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

