#' @export orphanMetabolites
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Identify the orphan metabolites of a set of stoichiometric reactions
#' @description This function identifies the orphan metabolites (metabolites not produced or not consumed in any other reaction or just involved in one reaction) for a set of stoichometric reactions.
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
#' @param actingAs A text string that specifies the type of metabolite to be returned; only \code{'reactant'} and \code{'product'} are supported.
#' @param byCompartment A boolean value \code{'TRUE'} or \code{'FALSE'} to indicate if orphan reactants should be reported by compartment
#' @return If \code{byCompartment == FALSE}, a vector with orphan reactants is returned, in opposite case a list is returned.
#' If \code{actingAs == 'reactant'}, metabolites not produced in any other reaction or just are involved in one reaction are returned.
#' If \code{actingAs == 'products'}, metabolites not consumed in any other reaction or just are involved in one reaction are returned.
#' @examples 
#' # Loading a set of stoichiometric reactions
#' glycolysis <- read.csv(system.file("extdata/glycolysisModel.csv",package = "minval"), sep="\t")
#' 
#' # Identify orphan metabolites
#' orphanMetabolites(reactionList = glycolysis$REACTION)
#' 
#' # Identify orphan reactants
#' orphanReactants(reactionList = glycolysis$REACTION)
#' 
#' # Identify orphan products
#' orphanProducts(reactionList = glycolysis$REACTION)
#' 
#' # Identify orphan metabolites by compartment
#' orphanMetabolites(reactionList = glycolysis$REACTION, byCompartment = TRUE)

orphanMetabolites <-
  function(reactionList,
           actingAs = NULL,
           byCompartment = FALSE) {
    # Convert to a vector
    reactionList <- as.vector(x = reactionList)
    # Remove reaction with invalid syntax
    reactionList <-
      reactionList[validateSyntax(reactionList = reactionList)]
    # Extract all reactants and products
    reactant <- unique(unlist(reactants(reactionList = reactionList)))
    product <- unique(unlist(products(reactionList = reactionList)))
    # Identify metabolites with frequency equal to 1
    met <-
      table(metabolites(reactionList = reactionList, uniques = FALSE))
    met <- names(met)[met == 1]
    # Identify orphans
    if (is.null(actingAs)) {
      orphan <-
        unique(c(reactant[!reactant %in% product], product[!product %in% reactant], met))
    } else if (actingAs == "reactant") {
      orphan <-
        unique(c(reactant[!reactant %in% product], met[met %in% reactant]))
    } else if (actingAs == "product") {
      orphan <-
        unique(c(product[!product %in% reactant], met[met %in% product]))
    } else {
      stop("actingAs should be one of: 'reactant' or 'product'")
    }
    # Return
    if (byCompartment == TRUE) {
      orphan <-
        sapply(compartments(orphan), function(compartment) {
          orphan[grep(paste0("\\[", compartment, "\\]"), orphan)]
        }, simplify = FALSE)
    }
    if (length(orphan) == 0){
      return (NA)
    } else {
      return(orphan)
    }
  }

#' @describeIn orphanMetabolites Identify the orphan reactants of a set of stoichometric reactions
#' @export orphanReactants
orphanReactants <- function(reactionList, byCompartment = FALSE) {
  orphanMetabolites(reactionList, actingAs = "reactant", byCompartment)
}

#' @describeIn orphanMetabolites Identify the orphan products of a set of stoichometric reactions
#' @export orphanProducts
orphanProducts <- function(reactionList, byCompartment = FALSE) {
  orphanMetabolites(reactionList, actingAs = "product", byCompartment)
}