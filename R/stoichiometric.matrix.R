stoichiometric.matrix <- function(reactionList){
  reactionList <- as.vector(reactionList)
  mets<- metabolites(reactionList)
  s <- matrix(0,nrow = length(reactionList),ncol=length(mets),dimnames = list(reactions=paste0("R",formatC(1:length(reactionList),digits = (nchar(length(reactionList))-1),flag = 0)),metabolites=mets))
  for (reaction in 1:length(reactionList)){
    r_met <- .get.reactant(reactionList[reaction])
    r_coe <- .coefficients(r_met)
    p_met <- .get.product(reactionList[reaction])
    p_coe <- .coefficients(p_met)
    r_met <- metabolites(r_met)
    p_met <- metabolites(p_met)
    s[reaction,r_met[!is.na(r_met)]] <- -1*(r_coe[!is.na(r_met)])
    s[reaction,p_met[!is.na(p_met)]] <- p_coe[!is.na(p_met)]
    s[reaction,r_met%in%p_met]<-0  
    
  }
  return(t(s))
}

