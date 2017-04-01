#' @export characterizeReactions
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Characterize stoichiometric reactions by compartments and reaction type
#' @description For a given set of stoichiometric reactions, this function: \itemize{
#' \item Counts the number of reactions, 
#' \item Computes the relative frequency of each reaction type (transport, exchange and compartmentalized), 
#' \item Computes the relative frequency of reactions by compartment, 
#' \item Counts the number of unique metabolites, 
#' \item Computes the relative frequency of metabolites by compartment.
#' }
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
#' @param rawOutput A boolean value \code{'TRUE'} or \code{'FALSE'} if computed values should be returned instead of raw data
#' @examples
#' # Loading a set of stoichiometric reactions
#' glycolysis <- read.csv(system.file("extdata/glycolysisModel.csv",package = "minval"), sep='\t')
#'
#' # Characterizing the reactions
#' characterizeReactions(reactionList = glycolysis$REACTION)
characterizeReactions <- function (reactionList, rawOutput = FALSE) {
  reactionList <-
    as.vector(reactionList)[validateSyntax(reactionList = reactionList)]
  model <- list ()
  model$nReactions <- length(unique(reactionList))
  model$rType <- reactionType(reactionList = reactionList)
  model$cReaction <-
    unlist(sapply(reactionList[model$rType == "Compartmentalized reaction"], function(reaction) {
      compartments(reactionList = reaction, uniques = TRUE)
    }, USE.NAMES = FALSE))
  model$rCompartments <- compartments(reactionList = reactionList)
  model$nMetabolites <-
    metabolites(reactionList = reactionList, uniques = TRUE)
  model$cMetabolites <-
    compartments(reactionList = model$nMetabolites, uniques = FALSE)
  model$nMetabolites <- length(model$nMetabolites)
  if (rawOutput == TRUE) {
    return(model)
  } else {
    summary <- list ()
    summary$nReactions <- model$nReactions
    summary$rType <- (table(model$rType) / model$nReactions) * 100
    summary$cReaction <-
      (table(model$cReaction) / model$nReactions) * 100
    summary$nMetabolites <- model$nMetabolites
    summary$cMetabolites <-
      (table(model$cMetabolites) / model$nMetabolites) * 100
    return(summary)
  }
}