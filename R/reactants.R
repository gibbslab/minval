#' @export reactants
#' @title Identify the reactants of a stoichometric reaction
#' @author Daniel Camilo Osorio <dcosorioh@tamu.edu>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @description This function identifies the reactants for a set of stoichometric reactions.
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
#' @return A vector with the identified reactants in the reaction, or a list with the identified reactats in each reaction if a set of stoichiometric reactions was given.

reactants <- function(reactionList) {
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Remove reaction with invalid syntax
  reactionList <- reactionList[validateSyntax(reactionList)]
  # Extract reactants for irreversible reactions
  reaction <-
    unlist(lapply(strsplit(reactionList, "[[:space:]]+=>[[:space:]]*"), function(reactionList) {
      reactionList[[1]]
    }))
  # Extract metabolites for reversible reactions
  reaction <- strsplit(reaction, "[[:space:]]+<=>[[:space:]]*")
  # Split independient reactants
  reaction <-
    lapply(reaction, function(reaction) {
      strsplit(reaction, "[[:space:]]+\\+[[:space:]]+")
    })
  # Remove spaces and report uniques
  reaction <-
    lapply(reaction, function(reaction) {
      unique(removeSpaces(unlist(reaction)))
    })
  # Use a regex to extract stoichiometric coefficients and separate the metabolite name
  reactants <-
    lapply(reaction, function(reaction) {
      unique(removeCoefficients(reaction))
    })
  if (length(reactants) == 0) {
    return(NA)
  } else if (length(reactants) == 1) {
    return(unlist(reactants))
  } else {
    return(reactants)
  }
}
