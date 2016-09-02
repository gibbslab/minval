#' @export isValidSyntax
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Evaluate if a stoichiometric reaction has a valid syntax
isValidSyntax <- function(reactionList){
  valid.syntax <- NULL
  # Coefficient validation
  valid.syntax <- c(valid.syntax,grepl("([[:digit:]][[:blank:]][[:digit:]][[:blank:]])+",reactionList))
  valid.syntax <- c(valid.syntax,grepl("(\\([[:digit:]]+\\)[[:blank:]]+)",reactionList))
  # Directionality validation
  valid.syntax <- c(valid.syntax, (!grepl("[[:blank:]]+<?=>[[:blank:]]*",reactionList)))
  valid.syntax <- c(valid.syntax,grepl("[[:blank:]](<?)-?-(>?)[[:blank:]]",reactionList))
  # Metabolite names validation
  valid.syntax <- c(valid.syntax,grepl("[[:alnum:]]+\\+[[:alnum:]]+",reactionList))
  # Blank spaces validation
  valid.syntax <- c(valid.syntax, (!grepl("[[:blank:]]",reactionList)))
  # Validating metabolite name
  valid.syntax <- c(valid.syntax, grepl("[[:blank:]]\\-[[:alnum:]]",reactionList))
  # Warnings!
  valid.syntax <- matrix(valid.syntax,ncol = 7)
  sapply(which(valid.syntax[,1]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid coefficients. Metabolites should have just one coefficient."),call. = FALSE)})
  sapply(which(valid.syntax[,2]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid coefficients. Coefficients should have not parentheses."),call. = FALSE)})
  sapply(which(valid.syntax[,3]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid directionality symbols. Arrow symbols should be between blank spaces."),call. = FALSE)})
  sapply(which(valid.syntax[,4]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid directionality symbols. Please use <=> or => instead of <-> or -> or -->."),call. = FALSE)})
  sapply(which(valid.syntax[,5]==TRUE),function(x){warning(paste0("Reaction ",x,": Metabolites names should be separated by a plus symbol between spaces."),call. = FALSE)})
  sapply(which(valid.syntax[,6]==TRUE),function(x){warning(paste0("Reaction ",x,": No blank spaces detected."),call. = FALSE)})
  sapply(which(valid.syntax[,7]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid metabolite name. Substituent should not be separated of the metabolite name."),call. = FALSE)})
  # Return
  valid.syntax <- sapply(seq_along(reactionList), function(reaction){identical(valid.syntax[reaction,],rep(FALSE,7))})
  return(valid.syntax)
}
