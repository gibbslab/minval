stoichiometric.matrix <- function(reactionList){
  reactionList <- as.vector(reactionList)
  mets<- metabolites(reactionList)
  s <- matrix(0,nrow = length(reactionList),ncol=length(mets),dimnames = list(reactions=paste0("R",formatC(1:length(reactionList),digits = (nchar(length(reactionList))-1),flag = 0)),metabolites=mets))
  for (reaction in 1:length(reactionList)){
    r_met <- .get.right(reactionList[reaction])
    r_coe <- .coefficients(r_met)
    p_met <- .get.left(reactionList[reaction])
    p_coe <- .coefficients(p_met)
    r_met <- metabolites(r_met)
    p_met <- metabolites(p_met)
    s[reaction,r_met] <- -1*(r_coe)
    s[reaction,p_met] <- p_coe
    s[reaction,r_met%in%p_met]<-0
  }
  return(t(s))
}

