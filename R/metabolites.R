#' @export metabolites
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Identify the list of metabolites for a set of stoichiometric reactions
#' @description This function identifies the list of metabolites for a set of stoichiometric reactions. If \code{'woCompartment'} is \code{'TRUE'} compartment label is removed. If \code{'uniques'} is \code{'TRUE'}, list of uniques is returned.
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
#' @param woCompartment A boolean value \code{'TRUE'} or \code{'FALSE'} to indicate if compartment label should be removed
#' @param uniques A boolean value \code{'TRUE'} or \code{'FALSE'} to indicate if uniques must be returned
#' @return A list of metabolites for a set of stoichiometric reactions
#' @examples 
#' # Extract metabolites of a stoichiometric reaction
#' metabolites(reactionList = "ADP[c] + Phosphoenolpyruvate[c] => ATP[c] + Pyruvate[c]")
#' 
#' # Loading a set of stoichiometric reactions
#' glycolysis <- read.csv(system.file("extdata/glycolysisModel.csv",package = "minval"), sep='\t')
#' 
#' # Extract unique metabolites
#' metabolites(reactionList = glycolysis$REACTION)
#' 
#' #' # Extract unique metabolites without compartments
#' metabolites(reactionList = glycolysis$REACTION, woCompartment = TRUE)
#' 
#' # Extract all metabolites
#' metabolites(reactionList = glycolysis$REACTION, uniques = FALSE)

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
