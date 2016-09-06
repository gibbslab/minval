#' @export compartments
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title  Extract the list of unique compartments for the metabolites of a set of stoichiometric reactions.
#' @description For a given set of stoichiometric reactions, this function identifies the compartments 
#' associated to each involved metabolite and return a vector with the list of unique compartments identified.
#' @param reactionList A set of stoichiometric reaction with the following format: \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"} Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}". 
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#' @return A vector with the list of of unique compartments identified for the metabolites of a set of stoichiometric reactions.
#' @examples
#' # Using individual reactions
#' compartments("H2O[c] + Urea-1-carboxylate[c] => 2 CO2[c] + 2 NH3[m]")
#' compartments("L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]")
#' 
#' # From a data.frame
#' glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
#' compartments(glycolysis$REACTION)
#' @keywords Extract Unique Compartments Metabolic Reconstruction

compartments <- function(reactionList){
  # Extract metabolites of a set of stoichiometric reactions
  metabolites <- metabolites(reactionList)
  # Extract compartments using a regular expression 
  compartments <- unique(unlist(regmatches(metabolites, gregexpr("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$", metabolites))))
  # Remove brackets
  compartments <- gsub("\\[|\\]","",compartments)
  # Return compartments
  return(compartments)
}
