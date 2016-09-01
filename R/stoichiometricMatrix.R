#' @export stoichiometricMatrix
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

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

