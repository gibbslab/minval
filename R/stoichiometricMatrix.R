#' @export stoichiometricMatrix
#' @title Build the stoichiometric matrix for a set of stoichiometric reactions
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

stoichiometricMatrix <- function(reactionList){
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Remove reaction with invalid syntax
  reactionList <- reactionList[validateSyntax(reactionList)]
  # Extract metabolites
  mets<- metabolites(reactionList,uniques = TRUE)
  # Create matrix
  s <- matrix(0,nrow = length(reactionList),ncol=length(mets),dimnames = list(reactions=paste0("R",formatC(1:length(reactionList),digits = (nchar(length(reactionList))-1),flag = 0)),metabolites=mets))
  # Fill
  for (reaction in seq_along(reactionList)){
    r_met <- unlist(getLeft(reactionList[reaction]))
    r_coe <- coefficients(r_met)
    p_met <- unlist(getRight(reactionList[reaction]))
    p_coe <- coefficients(p_met)
    r_met <- metabolites(r_met)
    p_met <- metabolites(p_met)
    if(any(r_met%in%p_met|p_met%in%r_met)){
      for(balanced in unique(c(r_met[r_met%in%p_met],p_met[p_met%in%r_met]))){
        if(r_coe[which(r_met == balanced)] > p_coe[which(p_met == balanced)]) {
          r_coe[which(r_met == balanced)] <- r_coe[which(r_met == balanced)] - p_coe[which(p_met == balanced)]
          p_coe <- p_coe[-which(p_met == balanced)]
          p_met <- p_met[-which(p_met == balanced)]
        } else if (r_coe[which(r_met == balanced)] < p_coe[which(p_met == balanced)]){
          p_coe[which(p_met == balanced)] <- p_coe[which(p_met == balanced)] - r_coe[which(r_met == balanced)]
          r_coe <- r_coe[-which(r_met == balanced)]
          r_met <- r_met[-which(r_met == balanced)]
        } else if (r_coe[which(r_met == balanced)] == p_coe[which(p_met == balanced)]){
          r_coe <- r_coe[-which(r_met == balanced)]
          p_coe <- p_coe[-which(p_met == balanced)]
          r_met <- r_met[-which(r_met == balanced)]
          p_met <- p_met[-which(p_met == balanced)]
        }
      }
    }
    s[reaction,r_met[!is.na(r_met)]] <- -1*(r_coe[!is.na(r_met)])
    s[reaction,p_met[!is.na(p_met)]] <- p_coe[!is.na(p_met)]
    s[reaction,r_met%in%p_met]<-0  
  }
  # Return
  return(t(s))
}
