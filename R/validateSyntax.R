#' @export validateSyntax
#' @author Daniel Camilo Osorio <dcosorioh@tamu.edu>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Evaluate if a stoichiometric reaction has a valid syntax
#' @description  For a set of given stoichiometric reactions, this function makes the following syntactic evaluations for each reaction: \itemize{
#' \item Evaluates if the reaction contain more than one coefficient by metabolite
#' \item Evaluates if the reaction contain metabolite coefficients between parenthesis
#' \item Evaluates if the reaction contain arrow symbol between spaces
#' \item Evaluates if the reaction contain not allowed arrow symbols
#' \item Evaluates if the reaction contain metabolites name separated by a plus symbol between spaces
#' \item Evaluates if the reaction contain substituents separated of the metabolite names
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
#' @return  A boolean value \code{'TRUE'} if reaction has a valid syntax.
#' @examples
#' # Evaluate the syntaxis for a single reaction
#' validateSyntax(reactionList = "ADP[c] + Phosphoenolpyruvate[c] => ATP[c] + Pyruvate[c]")
#'
#' # Loading a set of stoichiometric reactions
#' glycolysis <- read.csv(system.file("extdata/glycolysisModel.csv",package = "minval"), sep='\t')
#'
#' # Evaluating the syntaxis for a set of stoichiometric reactions
#' validateSyntax(reactionList = glycolysis$REACTION)

validateSyntax <- function(reactionList) {
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Create an empty vector
  valid.syntax <- NULL
  # Correct arrow symbol
  valid.syntax <-
    c(valid.syntax, (!(
      grepl("[[:space:]]+<=>[[:space:]]*", reactionList) |
        grepl("[[:space:]]+=>[[:space:]]*", reactionList)
    )))
  # Wrong arrow symbols
  valid.syntax <-
    c(valid.syntax, (
      grepl("[[:space:]]+<\\-[[:space:]]*", reactionList) |
        grepl("[[:space:]]+\\->[[:space:]]*", reactionList)
    ))
  valid.syntax <-
    c(valid.syntax, (
      grepl("[[:space:]]+<\\-\\-[[:space:]]*", reactionList) |
        grepl("[[:space:]]+\\-\\->[[:space:]]*", reactionList)
    ))
  valid.syntax <-
    c(valid.syntax, (
      grepl("[[:space:]]+<==>[[:space:]]*", reactionList) |
        grepl("[[:space:]]+==>[[:space:]]*", reactionList)
    ))
  # Wrong coefficients
  valid.syntax <-
    c(valid.syntax, (
      grepl(pattern = "([[:digit:]][[:space:]]+[[:digit:]][[:space:]]+)+", x = reactionList)
    ))
  # Coefficients surrounded by parentheses
  valid.syntax <-
    c(valid.syntax, (
      grepl(pattern = "(\\([[:digit:]]+\\)[[:space:]]+)", x = reactionList)
    ))
  # Correct metabolites distribution
  valid.syntax <-
    c(valid.syntax, (
      !grepl(pattern = "[[:graph:]]+[[:space:]]+<?=>[[:print:]]*", x = reactionList)
    ))
  # Metabolite name
  valid.syntax <-
    c(valid.syntax,
      grepl(pattern = "[[:space:]]+\\-[[:alnum:]]", x = reactionList))
  # Metabolite composition
  metabolitesList <-
    lapply(reactionList, function(reaction) {
      metabolites(reaction)
    })
  nMetabolites <- unlist(lapply(metabolitesList,function(reaction){length(reaction)}))
  valid.syntax <- c(valid.syntax, (nMetabolites < 1))
  valid.syntax <-
    c(valid.syntax, unlist(lapply(reactionList, function(reaction) {
      if (length(metabolites(reaction)) > 0) {
        left <- metabolites(unlist(getLeft(reaction)))
        left <- length(na.omit(left))
        right <- metabolites(unlist(getRight(reaction)))
        right <- length(na.omit(right))
        return(!(left >= 1 & right > 0 | left == 1 & right == 0))
      } else {
        return(FALSE)
      }
    })))
  valid.syntax <- c(
    valid.syntax, 
    unlist(lapply(metabolitesList,function(reaction){any(grepl("\\->",reaction) | grepl("<\\-",reaction))}))
    )
  
  # Warnings!
  valid.syntax <- matrix(valid.syntax, ncol = 11)
  sapply(which(valid.syntax[, 1] == TRUE), function(x) {
    warning(paste0("Reaction ", x, ": Not a valid arrow symbol was detected."),
            call. = FALSE)
  })
  sapply(which(valid.syntax[, 2] == TRUE), function(x) {
    warning(paste0("Reaction ", x, ": Not valid arrow symbol was detected."),
            call. = FALSE)
  })
  sapply(which(valid.syntax[, 3] == TRUE), function(x) {
    warning(paste0("Reaction ", x, ": Not valid arrow symbol was detected."),
            call. = FALSE)
  })
  sapply(which(valid.syntax[, 4] == TRUE), function(x) {
    warning(paste0("Reaction ", x, ": Not valid arrow symbol was detected."),
            call. = FALSE)
  })
  sapply(which(valid.syntax[, 5] == TRUE), function(x) {
    warning(
      paste0(
        "Reaction ",
        x,
        ": More than one stoichiometric coefficient by metabolite was detected."
      ),
      call. = FALSE
    )
  })
  sapply(which(valid.syntax[, 6] == TRUE), function(x) {
    warning(
      paste0(
        "Reaction ",
        x,
        ": Stoichiometric coefficients surrounded by parentheses were detected."
      ),
      call. = FALSE
    )
  })
  sapply(which(valid.syntax[, 7] == TRUE), function(x) {
    warning(
      paste0(
        "Reaction ",
        x,
        ": Metabolites not properly arranged in the biochemical reaction"
      ),
      call. = FALSE
    )
  })
  sapply(which(valid.syntax[, 8] == TRUE), function(x) {
    warning(
      paste0(
        "Reaction ",
        x,
        ": Substituents separated of the metabolite names were found"
      ),
      call. = FALSE
    )
  })
  sapply(which(valid.syntax[, 9] == TRUE), function(x) {
    warning(
      paste0(
        "Reaction ",
        x,
        ": Not metabolites were detected in biochemical reaction"
      ),
      call. = FALSE
    )
  })
  sapply(which(valid.syntax[, 10] == TRUE), function(x) {
    warning(
      paste0(
        "Reaction ",
        x,
        ": Metabolites not properly arranged in the biochemical reaction"
      ),
      call. = FALSE
    )
  })
  sapply(which(valid.syntax[, 11] == TRUE), function(x) {
    warning(
      paste0(
        "Reaction ",
        x,
        ": Metabolites must not contain arrows"
      ),
      call. = FALSE
    )
  })
  # Return
  if (any(valid.syntax == TRUE)) {
    message(
      "Please check that:\n* Arrows symbols are given in the form '=>' or '<=>'\n* Inverse arrow symbols '<=' or other types such as: '-->', '<==>' or '->' are not present\n* Arrows and plus signs are surrounded by a space character\n* Stoichiometric coefficients are surrounded by spaces and not by parentheses\n* Exchange reactions have only one metabolite before arrow symbol\n* Compartments are given between square brackets (metabolite[compartment]) joined at the end of metabolite name"
    )
  }
  return(sapply(seq_along(reactionList), function(reaction) {
    all(valid.syntax[reaction, ] == FALSE)
  }))
}
