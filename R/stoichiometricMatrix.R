#' @export stoichiometricMatrix
#' @title Return the stoichiometric matrix for a set of stoichiometric reactions
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @description  A set of stoichiometric reactions are often represented in a more compact form called the stoichiometry matrix. If a reaction network has n reactions and m participating metabolites then the stoichiometry matrix will have correspondingly m rows and n columns.
#' @param reactionList A set of stoichiometric reaction with the following format: 
#' 
#' \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"} 
#' 
#' Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}". 
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#' @return The stoichiometric matrix for a given set of stoichiometric reactions
#' @examples 
#' # Loading data
#' glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
#' 
#' # Removing stoichiometric reactions without valid syntax
#' glycolysis <- mapReactions(
#'                            reactionList = isValidSyntax(glycolysis$REACTION),
#'                            referenceData = glycolysis,
#'                            by = "bool"
#'                            )
#' # Building the Stoichiometric-Matrix 
#' stoichiometricMatrix(glycolysis$REACTION)
#' 
#' @keywords Stoichiometric Matrix Reactions Metabolic Reconstruction
stoichiometricMatrix <- function(reactionList){
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Remove reaction with invalid syntax
  reactionList <- reactionList[isValidSyntax(reactionList)]
  # Extract metabolites
  mets<- metabolites(reactionList)
  # Create matrix
  s <- matrix(0,nrow = length(reactionList),ncol=length(mets),dimnames = list(reactions=paste0("R",formatC(1:length(reactionList),digits = (nchar(length(reactionList))-1),flag = 0)),metabolites=mets))
  # Fill
  for (reaction in seq_along(reactionList)){
    r_met <- .get.left(reactionList[reaction])
    r_coe <- .coefficients(r_met)
    p_met <- .get.right(reactionList[reaction])
    p_coe <- .coefficients(p_met)
    r_met <- metabolites(r_met)
    p_met <- metabolites(p_met)
    s[reaction,r_met[!is.na(r_met)]] <- -1*(r_coe[!is.na(r_met)])
    s[reaction,p_met[!is.na(p_met)]] <- p_coe[!is.na(p_met)]
    s[reaction,r_met%in%p_met]<-0  
  }
  # Return
  return(t(s))
}

