# is.validsyntax
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

is.validsyntax <- function(reaction){
  valid.syntax <- NULL
  # Coefficient validation
  valid.syntax <- c(valid.syntax,grepl("([[:digit:]][[:blank:]][[:digit:]][[:blank:]])+",reaction))
  # Directionality validation
  valid.syntax <- c(valid.syntax, (!grepl("[[:blank:]]<*=>*[[:blank:]]",reaction)))
  valid.syntax <- c(valid.syntax,grepl("[[:blank:]](<?)-(>?)[[:blank:]]",reaction))
  # Metabolite names validation
  valid.syntax <- c(valid.syntax,grepl("[[:alnum:]]+\\+[[:alnum:]]+", reaction))
  # Blank spaces validation
  valid.syntax <- c(valid.syntax, (!grepl("[[:blank:]]",reaction)))
  #
  valid.syntax <- c(valid.syntax, grepl("[[:blank:]]\\-[[:alnum:]]",reaction))
  # Warnings!
  valid.syntax <- matrix(valid.syntax,ncol = 6)
  sapply(which(valid.syntax[,1]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid coefficients. Metabolites should have just one coefficient."),call. = FALSE)})
  sapply(which(valid.syntax[,2]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid directionality symbols. Arrow symbols should be between blank spaces."),call. = FALSE)})
  sapply(which(valid.syntax[,3]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid directionality symbols. Please use <=> or => instead of <-> or ->."),call. = FALSE)})
  sapply(which(valid.syntax[,4]==TRUE),function(x){warning(paste0("Reaction ",x,": Metabolites names should be separated by a plus symbol between spaces."),call. = FALSE)})
  sapply(which(valid.syntax[,5]==TRUE),function(x){warning(paste0("Reaction ",x,": No blank spaces detected."),call. = FALSE)})
  sapply(which(valid.syntax[,5]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid metabolite name. Substituent should not be separated of the metabolite name."),call. = FALSE)})
  # Return!
  valid.syntax <- sapply(1:length(reaction), function(reaction){identical(valid.syntax[reaction,],rep(FALSE,6))})
  return(valid.syntax)
}
