#' @export addExchangeReactions
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
addExchangeReactions <- function(reconstruction){
  reconstruction <- as.data.frame.array(reconstruction)
  if(all(names(reconstruction)%in%c("ID","DESCRIPTION","REACTION","GPR","LOWER.BOUND","UPPER.BOUND","OBJECTIVE"))){
    orphans <- table(c(orphanReactants(reconstruction[["REACTION"]]),orphanProducts(reconstruction[["REACTION"]])))
    new <- cbind(sapply(names(orphans), function(metabolite){paste0("EX_",metabolites(metabolite,woCompartment = TRUE),"(",compartments(metabolite),")")},USE.NAMES = FALSE),
                 sapply(names(orphans), function(metabolite){paste0("Exchange of ", metabolites(metabolite,woCompartment = TRUE))},USE.NAMES = FALSE),
                 sapply(names(orphans), function(metabolite){paste0(metabolite," <=>")},USE.NAMES = FALSE),
                 sapply(names(orphans), function(metabolite){""},USE.NAMES = FALSE),
                 ifelse(names(orphans)%in%orphanReactants(reconstruction[["REACTION"]]),-1000,0),
                 ifelse(names(orphans)%in%orphanProducts(reconstruction[["REACTION"]]),1000,0),
                 sapply(names(orphans), function(metabolite){0},USE.NAMES = FALSE))
    colnames(new) <- colnames(reconstruction)
    new <- as.data.frame(rbind(reconstruction,new),row.names=NULL)
    return(new)
  } else {
    return(FALSE)
  }
}