#' @export compartments
#' @author Daniel Camilo Osorio <dcosorioh@tamu.edu>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title  Extract the compartments associated to metabolites of a set of stoichiometric reactions.
#' @description For a given set of stoichiometric reactions, this function identifies the compartments
#' associated to each involved metabolite and return a vector with the list of unique compartments identified.
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
#' @param uniques A boolean value \code{'TRUE'} or \code{'FALSE'} if uniques must be returned
#' @return A vector with the list of compartments identified for the metabolites of a set of stoichiometric reactions.
#' @examples 
#' # Loading a set of stoichiometric reactions
#' glycolysis <- read.csv(system.file("extdata/glycolysisModel.csv",package = "minval"), sep='\t')
#' 
#' # Extract unique compartments
#' compartments(reactionList = glycolysis$REACTION)
#' 
#' # Extract all compartments
#' compartments(reactionList = glycolysis$REACTION, unique = FALSE)
#' 
#' # Extract compartments of metabolites
#' compartments(reactionList = "H2O[e]")

compartments <- function(reactionList, uniques = TRUE) {
  # Extract metabolites of a set of stoichiometric reactions
  if (uniques == TRUE) {
    metabolites <-
      metabolites(reactionList = reactionList, uniques = TRUE)
  } else {
    metabolites <-
      metabolites(reactionList = reactionList, uniques = FALSE)
  }
  # Extract compartments using a regular expression
  compartments <-
    unlist(regmatches(
      x = metabolites,
      m =  gregexpr(pattern = "\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$", text = metabolites)
    ))
  # Remove brackets
  compartments <- gsub("\\[|\\]", "", compartments)
  # NA return
  if (length(compartments) == 0) {
    compartments <- NA
  }
  # Return compartments
  if (uniques == TRUE) {
    compartments <- unique(compartments)
    return(compartments)
  } else {
    return(compartments)
  }
}
