#' @export stoichiometricMatrix
#' @title Build the stoichiometric matrix for a set of stoichiometric reactions
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @description  A set of stoichiometric reactions are often represented in a more compact form called the stoichiometry matrix.
#' If a metabolic network has n reactions and m participating metabolites then the stoichiometry matrix will have correspondingly m rows and n columns.
#' Values in stoichiometric matrix represent the metabolite coefficients in each reaction.
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
#' @return The stoichiometric matrix for a given set of stoichiometric reactions
#' @examples
#' # Loading a set of stoichiometric reactions
#' glycolysis <- read.csv(system.file("extdata/glycolysisModel.csv",package = "minval"), sep='\t')
#'
#' # Building the stoichiometric matrix
#' stoichiometricMatrix(reactionList = glycolysis$REACTION)

stoichiometricMatrix <- function(reactionList) {
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Remove reaction with invalid syntax
  reactionList <- reactionList[validateSyntax(reactionList)]
  # Extract metabolites
  mets <- metabolites(reactionList, uniques = TRUE)
  # Create matrix
  s <-
    matrix(
      0,
      nrow = length(reactionList),
      ncol = length(mets),
      dimnames = list(
        reactions = paste0("R", formatC(
          1:length(reactionList),
          digits = (nchar(length(reactionList)) - 1),
          flag = 0
        )),
        metabolites = mets
      )
    )
  # Fill
  for (reaction in seq_along(reactionList)) {
    r_met <- unlist(getLeft(reactionList[reaction]))
    r_coe <- coefficients(r_met)
    p_met <- unlist(getRight(reactionList[reaction]))
    p_coe <- coefficients(p_met)
    r_met <- metabolites(r_met)
    p_met <- metabolites(p_met)
    if (any(r_met %in% p_met | p_met %in% r_met)) {
      for (balanced in unique(c(r_met[r_met %in% p_met], p_met[p_met %in% r_met]))) {
        if (r_coe[which(r_met == balanced)] > p_coe[which(p_met == balanced)]) {
          r_coe[which(r_met == balanced)] <-
            r_coe[which(r_met == balanced)] - p_coe[which(p_met == balanced)]
          p_coe <- p_coe[-which(p_met == balanced)]
          p_met <- p_met[-which(p_met == balanced)]
        } else if (r_coe[which(r_met == balanced)] < p_coe[which(p_met == balanced)]) {
          p_coe[which(p_met == balanced)] <-
            p_coe[which(p_met == balanced)] - r_coe[which(r_met == balanced)]
          r_coe <- r_coe[-which(r_met == balanced)]
          r_met <- r_met[-which(r_met == balanced)]
        } else if (r_coe[which(r_met == balanced)] == p_coe[which(p_met == balanced)]) {
          r_coe <- r_coe[-which(r_met == balanced)]
          p_coe <- p_coe[-which(p_met == balanced)]
          r_met <- r_met[-which(r_met == balanced)]
          p_met <- p_met[-which(p_met == balanced)]
        }
      }
    }
    s[reaction, r_met[!is.na(r_met)]] <- -1 * (r_coe[!is.na(r_met)])
    s[reaction, p_met[!is.na(p_met)]] <- p_coe[!is.na(p_met)]
    s[reaction, r_met %in% p_met] <- 0
  }
  # Return
  return(t(s))
}
