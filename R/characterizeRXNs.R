#' @export characterizeRXNs
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
characterizeRXNs <- function (reactionList, rawOutput = FALSE){
  reactionList <- as.vector(reactionList)[validateSyntax(reactionList = reactionList)]
  model <- list ()
  model$nReactions <- length(unique(reactionList))
  model$rType <- reactionType(reactionList = reactionList)
  model$cReaction <- unlist(sapply(reactionList[model$rType == "Compartmentalized reaction"], function(reaction){compartments(reactionList = reaction, uniques = TRUE)},USE.NAMES = FALSE))
  model$rCompartments <- compartments(reactionList = reactionList)
  model$nMetabolites <- metabolites(reactionList = reactionList, uniques = TRUE)
  model$cMetabolites <- compartments(reactionList = model$nMetabolites, uniques = FALSE)
  model$nMetabolites <- length(model$nMetabolites)
  if (rawOutput == TRUE){
    return(model)
  } else {
    summary <- list ()
    summary$nReactions <- model$nReactions
    summary$rType <- (table(model$rType) / model$nReactions) * 100
    summary$cReaction <- (table(model$cReaction)/ model$nReactions) * 100
    summary$nMetabolites <- model$nMetabolites
    summary$cMetabolites <- (table(model$cMetabolites)/ model$nMetabolites) * 100
    return(summary)
  }
}